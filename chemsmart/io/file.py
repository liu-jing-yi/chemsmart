import logging
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

import numpy as np

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.io import (
    _reposition_rings_and_metal,
    attach_eta_bonds_for_arene_rings,
    attach_eta_bonds_for_cp_rings,
)
from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


def _find_element_parent(root, target):
    """Return the parent element of *target* inside *root*, or None."""
    for parent in root.iter():
        for child in parent:
            if child is target:
                return parent
    return None


def _cdxml_atom_color_counts(atoms):
    """Return colour-table counts for node colours and implicit-H spans."""
    colors = []
    for atom in atoms:
        colors.append(atom["color"])
        if atom["implicit_h_color"] is not None:
            colors.append(atom["implicit_h_color"])
    return Counter(colors)


_PKB_COLORED_HYDROGEN_ERROR = (
    "Uniquely coloured atom is hydrogen. With --pkb, colour the basic "
    "heavy atom to protonate, or drop --pkb and run pKa."
)


class SDFFile(FileMixin):
    """
    SDF file object.
    """

    def __init__(self, filename):
        self.filename = filename

    @property
    def molecule(self):
        return self.get_molecule()

    def get_molecule(self):
        list_of_symbols = []
        cart_coords = []
        # sdf line pattern containing coordinates and element type
        from chemsmart.utils.repattern import sdf_pattern

        for line in self.contents:
            match = re.match(sdf_pattern, line)
            if match:
                x = float(match.group(1))
                y = float(match.group(2))
                z = float(match.group(3))
                atom_type = str(match.group(4))
                list_of_symbols.append(atom_type)
                cart_coords.append((x, y, z))

        cart_coords = np.array(cart_coords)

        if len(list_of_symbols) == 0 or len(cart_coords) == 0:
            raise ValueError("No coordinates found in the SDF file!")

        return Molecule.from_symbols_and_positions_and_pbc_conditions(
            list_of_symbols=list_of_symbols, positions=cart_coords
        )


class CDXFile(FileMixin):
    """
    ChemDraw file object for reading .cdx and .cdxml files.

    Supports both binary (.cdx) and XML-based (.cdxml) ChemDraw formats.
    Uses RDKit for parsing and generates 3D coordinates using EmbedMolecule.

    Args:
        filename (str or pathlib.Path): Path to the ChemDraw file.
    """

    def __init__(self, filename):

        # Accept both str and pathlib.Path
        self.filename = (
            str(filename) if isinstance(filename, Path) else filename
        )

    @property
    def molecules(self):
        """
        Return all molecules from the ChemDraw file.
        """
        return self._parse_chemdraw_file()

    def get_molecules(self, index="-1", return_list=False):
        """Return molecule(s) for the requested 1-based ChemDraw index."""
        molecules = self.molecules

        if index == ":":
            selection = molecules
        else:
            parsed_index = (
                index
                if isinstance(index, slice)
                else string2index_1based(str(index))
            )
            if isinstance(parsed_index, slice):
                selection = molecules[parsed_index]
            else:
                selection = molecules[parsed_index]

        if return_list:
            return selection if isinstance(selection, list) else [selection]
        return selection

    def _parse_chemdraw_file(self):
        """
        Parse the ChemDraw file and return a list of Molecule objects.

        Uses RDKit to parse the file and generate 3D coordinates.
        Falls back to Open Babel (via obtain_mols_from_cdx_via_obabel) for .cdx
        if RDKit cannot read it.

        Returns:
            list[Molecule]: List of Molecule objects with 3D coordinates.

        Raises:
            ValueError: If no molecules can be read or no valid 3D structures
            can be generated from the ChemDraw file.
        """

        from rdkit import Chem

        suffix = Path(self.filename).suffix.lower()

        rdkit_mols = []
        try:
            # NOTE: RDKit's MolsFromCDXMLFile always supports CDXML.
            # CDX files are only supported if RDKit was built with ChemDraw CDX support.
            # Use sanitize=False to avoid kekulization errors during parsing
            # of organometallic complexes. We'll handle sanitization later.
            logger.debug(
                f"Generating rdkit mols from {self.filename} using RDKit"
            )
            if suffix == ".cdxml":
                # Preprocess CDXML to remove MultiAttachment phantom atoms
                # BEFORE RDKit reads the file.  ChemDraw represents η5/η6
                # hapticity with NodeType="MultiAttachment" nodes connected to
                # the metal via Display="Dash" bonds.  RDKit turns each such
                # node into a regular carbon atom, producing spurious CH₃
                # groups and masking real substituents (e.g. Ti–Me groups).
                # Removing them at the XML level keeps genuine ligands intact.
                logger.debug(
                    f"Preprocessing CDXML to remove MultiAttachment nodes: {self.filename}"
                )
                cleaned_cdxml = (
                    CDXFile._preprocess_cdxml_remove_multi_attachments(
                        self.filename
                    )
                )
                rdkit_mols = list(
                    Chem.MolsFromCDXML(
                        cleaned_cdxml, sanitize=False, removeHs=False
                    )
                )
            else:
                rdkit_mols = list(
                    Chem.MolsFromCDXMLFile(
                        self.filename, sanitize=False, removeHs=False
                    )
                )
        except Exception as e:
            logger.debug(
                f"RDKit MolsFromCDXMLFile failed for {self.filename}: {e}"
            )

        # Fallback for .cdx: use Open Babel helper if RDKit gave nothing.
        if not rdkit_mols and suffix == ".cdx":
            logger.debug(
                f"RDKit did not return any molecules for {self.filename}; "
                "falling back to Open Babel CDX reader."
            )
            try:
                from chemsmart.utils.io import obtain_mols_from_cdx_via_obabel

                logger.debug(
                    f"Generating rdkit mols from {self.filename} using obabel"
                )
                rdkit_mols = obtain_mols_from_cdx_via_obabel(self.filename)
            except Exception as e:
                logger.debug(
                    f"Open Babel CDX fallback failed for {self.filename}: {e}"
                )

        if not rdkit_mols:
            raise ValueError(
                f"No molecules could be read from ChemDraw file: {self.filename}"
            )

        # Update property cache for all molecules read from CDXML
        # This is necessary to avoid "Pre-condition Violation" errors
        # when checking atom properties like GetTotalNumHs()
        for mol in rdkit_mols:
            if mol is not None:
                mol.UpdatePropertyCache(strict=False)

        # Combine metal fragments with their aromatic ligands
        # ChemDraw sometimes draws metal complexes as separate fragments
        logger.debug("Combining metal fragments with ligands")
        rdkit_mols = self._combine_metal_and_ligand_fragments(rdkit_mols)

        molecules = []
        for rdkit_mol in rdkit_mols:
            if rdkit_mol is None:
                continue

            # Process molecule (handle organometallic complexes, add H, generate 3D coords)
            try:
                logger.debug(f"Processing rdkit mol: {rdkit_mol}")
                rdkit_mol = self._process_cdx_molecule(rdkit_mol)
            except Exception as e:
                logger.warning(
                    f"Error processing molecule in {self.filename}: {e}. "
                    f"Skipping this molecule."
                )
                continue

            # Convert RDKit Mol to Molecule
            mol = Molecule.from_rdkit_mol(rdkit_mol)
            molecules.append(mol)

        if not molecules:
            raise ValueError(
                "No valid molecules with 3D coordinates could be generated "
                f"from ChemDraw file: {self.filename}"
            )

        return molecules

    def _process_cdx_molecule(self, rdkit_mol):
        """
        Process a molecule from ChemDraw file, handling organometallic complexes.

        This method normalizes metal bonds and sanitizes the molecule appropriately
        for both organic and organometallic structures. It then adds hydrogens and
        generates 3D coordinates.

        Args:
            rdkit_mol (rdkit.Chem.Mol): RDKit molecule to process.

        Returns:
            rdkit.Chem.Mol: Processed molecule with hydrogens and 3D coordinates.

        Raises:
            Exception: If molecule processing fails.
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem

        from chemsmart.utils.io import normalize_metal_bonds, safe_sanitize
        from chemsmart.utils.periodictable import NON_METALS_AND_METALLOIDS

        # Check if molecule contains metals
        has_metals = any(
            atom.GetAtomicNum() not in NON_METALS_AND_METALLOIDS
            for atom in rdkit_mol.GetAtoms()
        )
        logger.debug(
            f"The molecule: {rdkit_mol} has {has_metals} metal atoms."
        )

        logger.debug(f"Building metal indices for {rdkit_mol}.")
        metal_idxs = {
            a.GetIdx()
            for a in rdkit_mol.GetAtoms()
            if a.GetAtomicNum() not in NON_METALS_AND_METALLOIDS
        }

        # Normalize metal bonds first (removes aromatic flags from metal bonds)
        logger.debug(f"Normalize metal bonds in {rdkit_mol}.")
        rdkit_mol = normalize_metal_bonds(rdkit_mol)

        if has_metals:
            # Add ONE bond from the metal to each isolated 5-membered Cp ring
            # and each isolated 6-membered arene ring.  This also dearomatizes
            # the Cp ring with alternating single/double bonds so every ring
            # carbon is sp2 with exactly one implicit H.  One bond per ring is
            # enough for ETKDG to embed the molecule; the metal is repositioned
            # after embedding.
            logger.debug(f"Attach η5 bonds for Cp rings in {rdkit_mol}.")
            rdkit_mol = attach_eta_bonds_for_cp_rings(rdkit_mol, metal_idxs)
            logger.debug(f"Attach η6 bonds for arene rings in {rdkit_mol}.")
            rdkit_mol = attach_eta_bonds_for_arene_rings(rdkit_mol, metal_idxs)

        # Sanitize with or without kekulization based on metal presence
        logger.debug(f"Sanitize metal bonds in {rdkit_mol}.")
        rdkit_mol = safe_sanitize(rdkit_mol, skip_kekulize=has_metals)

        # Add explicit hydrogens for proper structure
        logger.debug(f"Adding explicit hydrogens in {rdkit_mol}.")
        rdkit_mol = Chem.AddHs(rdkit_mol)

        # Generate 3D coordinates
        # Try to embed the molecule to get 3D coordinates
        logger.debug(f"Embed molecule in {rdkit_mol}.")
        result = AllChem.EmbedMolecule(rdkit_mol, randomSeed=42)
        if result == -1:
            # Embedding failed, try with random coordinates
            logger.debug(f"Failed to embed molecule in {rdkit_mol}.")
            result = AllChem.EmbedMolecule(
                rdkit_mol,
                useRandomCoords=True,
                randomSeed=42,
            )
            if result == -1:
                raise ValueError(
                    "Could not generate 3D coordinates for molecule"
                )

        # For organometallic complexes, ETKDG places the metal and all ring
        # atoms in the same plane (coplanar with the single anchor bond).
        # Reposition ring atoms as rigid bodies and then move the metal to
        # achieve proper η5/η6 sandwich/half-sandwich coordination geometry.
        if has_metals:
            logger.debug(f"Reposition rings and metal in {rdkit_mol}.")
            rdkit_mol = _reposition_rings_and_metal(rdkit_mol, metal_idxs)

        # Optimize the geometry (may fail for exotic atom types)
        try:
            logger.debug(f"Optimize molecule {rdkit_mol} using MMFF.")
            AllChem.MMFFOptimizeMolecule(rdkit_mol)
        except Exception as e:
            logger.debug(
                f"MMFF optimization failed for a molecule in {self.filename}: {e}"
            )

        return rdkit_mol

    def _combine_metal_and_ligand_fragments(self, rdkit_mols):
        """
        Combine metal fragments with their aromatic ligand fragments.

        ChemDraw sometimes draws organometallic complexes as separate fragments
        (e.g., a metal center with coordinated aromatic rings as separate fragments).
        This method detects such cases and combines them into single molecules.

        The heuristic used:
        - A small metal-containing fragment (< 10 atoms) followed by aromatic ring
          fragments (6-atom benzene rings) are combined.

        Args:
            rdkit_mols (list[rdkit.Chem.Mol]): List of RDKit molecules (fragments).

        Returns:
            list[rdkit.Chem.Mol]: List of combined molecules.
        """
        from rdkit import Chem

        from chemsmart.utils.periodictable import NON_METALS_AND_METALLOIDS

        def has_metal(mol):
            """Check if molecule contains metal atoms."""
            if mol is None:
                return False
            return any(
                atom.GetAtomicNum() not in NON_METALS_AND_METALLOIDS
                for atom in mol.GetAtoms()
            )

        def is_small_ligand_ring(mol):
            """Check if molecule is a small aromatic ring (e.g., benzene, Cp)."""
            if mol is None:
                return False
            # Check for aromatic 5-member carbon ring (Cp) or 6-member carbon ring (benzene)
            ri = mol.GetRingInfo()
            if ri.NumRings() != 1:
                return False
            for ring in ri.AtomRings():
                if len(ring) in (5, 6):
                    ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                    if all(a.GetSymbol() == "C" for a in ring_atoms):
                        return True
            return False

        combined_mols = []
        i = 0
        while i < len(rdkit_mols):
            mol = rdkit_mols[i]

            if mol is None:
                i += 1
                continue

            # Check if this is a small metal-containing fragment
            if has_metal(mol) and mol.GetNumAtoms() < 10:
                # Look ahead for aromatic ligand fragments
                ligands = []
                j = i + 1
                while j < len(rdkit_mols) and is_small_ligand_ring(
                    rdkit_mols[j]
                ):
                    ligands.append(rdkit_mols[j])
                    j += 1

                # If we found ligands, combine them with the metal fragment
                if ligands:
                    logger.debug(
                        f"Combining metal fragment {i} with {len(ligands)} ligand(s)"
                    )
                    # Remove degree-1 carbon stub atoms from the small metal
                    # fragment before combining.  ChemDraw sometimes draws
                    # η5/η6 complexes as a small metal fragment (M + 2 C stubs)
                    # plus separate ring fragments.  The C stubs are drawing
                    # artefacts representing the bond direction to the ring
                    # centroid; they are not real atoms and must be removed so
                    # that attach_eta_bonds_for_{cp,arene}_rings can later
                    # attach the metal directly to the ring atoms.
                    stub_idxs = [
                        nbr.GetIdx()
                        for a in mol.GetAtoms()
                        if a.GetAtomicNum() not in NON_METALS_AND_METALLOIDS
                        for nbr in a.GetNeighbors()
                        if nbr.GetAtomicNum() == 6 and nbr.GetDegree() == 1
                    ]
                    if stub_idxs:
                        rw = Chem.RWMol(mol)
                        for idx in sorted(set(stub_idxs), reverse=True):
                            rw.RemoveAtom(idx)
                        mol = rw.GetMol()
                        mol.UpdatePropertyCache(strict=False)

                    combined = Chem.CombineMols(mol, ligands[0])
                    for k in range(1, len(ligands)):
                        combined = Chem.CombineMols(combined, ligands[k])

                    # Update property cache after combining molecules
                    # This is critical to avoid "Pre-condition Violation" errors
                    combined.UpdatePropertyCache(strict=False)
                    combined_mols.append(combined)
                    i = j  # Skip the ligands we just combined
                else:
                    combined_mols.append(mol)
                    i += 1
            else:
                combined_mols.append(mol)
                i += 1

        return combined_mols

    @staticmethod
    def _preprocess_cdxml_remove_multi_attachments(filepath):
        """
        Read a CDXML file and return the XML content with all
        ``NodeType="MultiAttachment"`` atoms and their connecting bonds removed.

        ChemDraw uses ``MultiAttachment`` nodes to represent η5/η6 hapticity
        in structures like metallocenes.  These phantom atoms are connected to
        the metal via ``Display="Dash"`` bonds and carry no chemical meaning.
        When RDKit reads them as regular carbon atoms it produces spurious CH₃
        groups *and* masks real substituents bonded to the metal (e.g. the two
        Ti–Me groups in TiCp₂Me₂).

        Removing the phantom atoms at the XML level—before any RDKit
        processing—preserves every genuine ligand while leaving the ring atoms
        as isolated (disconnected) components inside the same CDXML fragment.
        ``attach_eta_bonds_for_cp_rings`` / ``attach_eta_bonds_for_arene_rings``
        then reconnect the metal to the rings in the correct η5/η6 fashion.

        Args:
            filepath (str): Path to the CDXML file.

        Returns:
            str: Modified CDXML content as a string (DOCTYPE declaration
            stripped, which is acceptable for ``Chem.MolsFromCDXML``).
        """
        import io
        import xml.etree.ElementTree as ET

        tree = ET.parse(filepath)
        root = tree.getroot()

        for frag in root.iter("fragment"):
            # Collect IDs of MultiAttachment atoms in this fragment.
            multi_ids: set[str] = set()
            nodes_to_remove = []
            for child in list(frag):
                if (
                    child.tag == "n"
                    and child.get("NodeType") == "MultiAttachment"
                ):
                    node_id = child.get("id")
                    if node_id:
                        multi_ids.add(node_id)
                    nodes_to_remove.append(child)

            if not multi_ids:
                continue

            # Collect bonds that reference any MultiAttachment atom.
            bonds_to_remove = [
                child
                for child in list(frag)
                if child.tag == "b"
                and (
                    child.get("B") in multi_ids or child.get("E") in multi_ids
                )
            ]

            for node in nodes_to_remove:
                frag.remove(node)
            for bond in bonds_to_remove:
                frag.remove(bond)

        buf = io.StringIO()
        tree.write(buf, encoding="unicode", xml_declaration=False)
        return buf.getvalue()

    # ------------------------------------------------------------------
    # CDXML atom-colour helpers
    # ------------------------------------------------------------------


class PKaCDXFile(CDXFile):
    """Specialized ``CDXFile`` subclass used in pKa workflows.

    This class adds utilities on top of :class:`CDXFile` for:

    * parsing ChemDraw CDXML atom colours and element labels,
    * distinguishing explicitly drawn hydrogens from implicit/functional-group
      protons based on label structure and colour differences, and
    * providing structured atom and proton information that can be consumed
      by higher-level pKa logic (for example, to build :class:`PKaMolecule`
      instances or similar pKa-aware molecule objects).

    The class does **not** perform the pKa calculation itself; instead it
    focuses on converting CDXML drawings into a representation where
    ionizable sites and their associated protons can be identified reliably.
    """

    def _resolve_proton_from_cdxml(self, color_code=None, mode="acid"):
        """Auto-detect the ionizable site from CDX/CDXML colour markup.

        Args:
            color_code: CDXML color-table index used for auto-detection.
            mode: ``"acid"`` (coloured hydrogen) or ``"base"`` (coloured
                heavy atom to protonate).

        Returns:
            tuple[int | None, list | None]:
                - Site index when a single molecule is resolved.
                - ``None`` for the index with a list of per-fragment molecules
                  when multiple molecules are detected.

        Raises:
            ValueError: If CDXML colour parsing fails.
        """
        if mode == "base":
            try:
                pairs = list(
                    self._iter_colored_sites(
                        color_code=color_code, mode="base"
                    )
                )
            except ValueError as exc:
                raise ValueError(
                    "Could not auto-detect basic atom from CDXML colour: "
                    f"{exc}\n"
                    "Use -pi/--proton-index to specify the atom explicitly."
                ) from exc
            if len(pairs) > 1:
                mols = []
                for mol, site in pairs:
                    mol.proton_index = site
                    mols.append(mol)
                logger.info(
                    "Detected %s molecules with per-fragment basic-atom "
                    "auto-detection in %s.",
                    len(mols),
                    self.filename,
                )
                return None, mols
            site = pairs[0][1]
            logger.info(
                "Detected basic atom index %s from CDXML colour in %s.",
                site,
                self.filename,
            )
            return site, None

        try:
            pka_mols = self.get_pka_molecules(
                color_code=color_code, index=":", return_list=True
            )
        except ValueError as exc:
            raise ValueError(
                "Could not auto-detect proton from CDXML colour: "
                f"{exc}\n"
                "Use -pi/--proton-index to specify the proton explicitly."
            ) from exc

        if len(pka_mols) > 1:
            logger.info(
                "Detected %s molecules with per-fragment proton "
                "auto-detection in %s.",
                len(pka_mols),
                self.filename,
            )
            return None, pka_mols

        proton_index = pka_mols[0].proton_index
        logger.info(
            "Detected proton index %s from CDXML colour in %s.",
            proton_index,
            self.filename,
        )
        return proton_index, None

    @staticmethod
    def resolve_reference_proton(
        reference, reference_proton_index, reference_color_code
    ):
        """Resolve reference proton index for CDX/CDXML files.

        Args:
            reference: Reference structure file path.
            reference_proton_index: 1-based proton index supplied by user.
            reference_color_code: CDXML color-table index used for detection.

        Returns:
            int | None: Resolved reference proton index, or None when reference
                is not provided.

        Raises:
            ValueError: If CDXML parsing fails or color-code usage is invalid.
        """
        if reference_proton_index is not None:
            return reference_proton_index
        if reference is None:
            return None

        reference = str(reference)
        if reference.endswith((".cdx", ".cdxml")):
            ref_cdx = PKaCDXFile(reference)
            ref_pka_mol = ref_cdx.get_pka_molecules(
                index="-1", color_code=reference_color_code
            )
            return ref_pka_mol.proton_index

        if reference_color_code is not None:
            raise ValueError(
                "-rcc/--reference-color-code can only be used when --reference "
                "is a .cdx/.cdxml file."
            )

        return None

    def parse_cdxml_element_colors(self):
        """Parse the CDXML file and return per-atom colour information.

        Walks every ``<fragment>`` element in document order and collects
        ``<n>`` (node / atom) elements.  Atoms without an ``Element``
        attribute are implicit carbons and are still tracked for correct
        index mapping.

        For each atom the method also inspects the ``<t>``/``<s>`` text
        spans that make up the displayed label.  When a label contains
        an "H" character rendered in a **different** colour from the
        heavy-atom symbol (e.g. ``<s color="0">O</s><s color="4">H</s>``
        for a phenol –OH), the colour of that "H" span is recorded in
        ``implicit_h_color``.  This enables detection of functional-group
        protons (–OH, –NH, …) that are *not* represented as separate
        ``<n>`` nodes but instead as implicit hydrogens on the heavy atom.

        Returns:
            list[dict]: One entry per CDXML atom, in document order, with
                keys:

                * ``cdxml_id`` – the ``id`` attribute of the ``<n>`` element
                * ``element`` – integer atomic number (6 for C when omitted)
                * ``color`` – integer colour-table index (``0`` when absent)
                * ``symbol`` – element symbol string (``"C"``, ``"H"``, …)
                * ``num_hydrogens`` – number of implicit hydrogens attached
                  (from the ``NumHydrogens`` CDXML attribute, or ``None``)
                * ``implicit_h_color`` – colour-table index of an "H" span
                  in the label that differs from the atom's own colour, or
                  ``None`` when no such span exists

        Raises:
            ValueError: If the file cannot be parsed as valid CDXML.

        Note:
            The atom order matches the order RDKit uses when reading the
            same CDXML through ``Chem.MolsFromCDXMLFile``, **before**
            ``Chem.AddHs`` adds implicit hydrogens.
        """
        from chemsmart.utils.periodictable import PeriodicTable

        pt = PeriodicTable()
        root = self._parse_cdxml_root()

        list_of_elements = []
        for fragment in self._iter_top_level_fragments(root):
            list_of_elements.extend(self._parse_fragment_nodes(fragment, pt))

        return list_of_elements

    def _parse_cdxml_root(self):
        """Parse CDXML and return the XML root element."""
        try:
            tree = ET.parse(self.filename)
        except ET.ParseError as exc:
            raise ValueError(
                f"Failed to parse CDXML file {self.filename}: {exc}"
            ) from exc
        return tree.getroot()

    def _iter_top_level_fragments(self, root):
        """Yield top-level fragment elements, excluding nested node fragments."""
        for fragment in root.iter("fragment"):
            parent = _find_element_parent(root, fragment)
            if parent is not None and parent.tag == "n":
                continue
            yield fragment

    def _parse_fragment_nodes(self, fragment, periodic_table):
        """Convert all valid CDXML nodes in one fragment into atom dicts."""
        fragment_atoms = []
        for node in fragment.findall("n"):
            if node.get("NodeType") == "ExternalConnectionPoint":
                continue

            cdxml_id = node.get("id")
            element_num = int(node.get("Element", "6"))
            color = int(node.get("color", "0"))
            num_h_attr = node.get("NumHydrogens")
            num_hydrogens = int(num_h_attr) if num_h_attr is not None else None

            spans = []
            for t_elem in node.iter("t"):
                for s_elem in t_elem.iter("s"):
                    s_text = (s_elem.text or "").strip()
                    s_color = int(s_elem.get("color", "0"))
                    if s_text:
                        spans.append((s_text, s_color))

            if color == 0 and spans:
                for s_text, s_color in spans:
                    if s_text.upper() not in ("H",):
                        color = s_color
                        break
                else:
                    color = spans[0][1]

            implicit_h_color = None
            for s_text, s_color in spans:
                if "H" in s_text and s_color != color:
                    implicit_h_color = s_color
                    break

            try:
                symbol = periodic_table.to_symbol(element_num)
            except Exception:
                symbol = "?"

            fragment_atoms.append(
                {
                    "cdxml_id": cdxml_id,
                    "element": element_num,
                    "color": color,
                    "symbol": symbol,
                    "num_hydrogens": num_hydrogens,
                    "implicit_h_color": implicit_h_color,
                }
            )

        return fragment_atoms

    def _rdkit_atom_idx_by_cdxml_id(self, rdkit_mol, cdxml_id):
        """Return the 0-based RDKit index for a CDXML node id."""
        if rdkit_mol is None:
            return None
        target_id = int(cdxml_id)
        for atom in rdkit_mol.GetAtoms():
            if atom.HasProp("_CDX_ATOM_ID"):
                if atom.GetIntProp("_CDX_ATOM_ID") == target_id:
                    return atom.GetIdx()
        return None

    def _nested_deprotonatable_cdxml_ids(self, node_cdxml_id):
        """Return nested atom ids with explicit/implicit H inside a fragment node."""
        root = self._parse_cdxml_root()
        node = None
        for candidate in root.iter("n"):
            if candidate.get("id") == str(node_cdxml_id):
                node = candidate
                break
        if node is None or node.get("NodeType") != "Fragment":
            return []

        inner_fragment = node.find("fragment")
        if inner_fragment is None:
            return []

        ids = []
        for child in inner_fragment.findall("n"):
            if child.get("NodeType") == "ExternalConnectionPoint":
                continue
            num_h_attr = child.get("NumHydrogens")
            if num_h_attr is None or int(num_h_attr) < 1:
                continue
            element = int(child.get("Element", "0"))
            if element in (1, 7, 8, 16):
                child_id = child.get("id")
                if child_id is not None:
                    ids.append(child_id)
        return ids

    def _rdkit_heavy_idx_for_implicit_h(
        self, rdkit_mol_h, atom, local_idx=None
    ):
        """Map a coloured implicit-H heavy atom to a post-AddHs RDKit index."""
        heavy_idx = self._rdkit_atom_idx_by_cdxml_id(
            rdkit_mol_h, atom["cdxml_id"]
        )
        if heavy_idx is None:
            for nested_id in self._nested_deprotonatable_cdxml_ids(
                atom["cdxml_id"]
            ):
                heavy_idx = self._rdkit_atom_idx_by_cdxml_id(
                    rdkit_mol_h, nested_id
                )
                if heavy_idx is not None:
                    break

        if heavy_idx is None and local_idx is not None:
            if 0 <= local_idx < rdkit_mol_h.GetNumAtoms():
                heavy_idx = local_idx
        return heavy_idx

    def _proton_index_from_rdkit_heavy(self, rdkit_mol_h, heavy_idx, atom):
        """Return the 1-based index of an H bonded to *heavy_idx* after AddHs."""
        heavy_atom = rdkit_mol_h.GetAtomWithIdx(heavy_idx)
        h_indices = [
            bond.GetOtherAtom(heavy_atom).GetIdx()
            for bond in heavy_atom.GetBonds()
            if bond.GetOtherAtom(heavy_atom).GetSymbol() == "H"
        ]
        if not h_indices:
            raise ValueError(
                f"No hydrogen bonded to {atom['symbol']} "
                f"(CDXML id={atom['cdxml_id']}) after AddHs."
            )

        proton_index = h_indices[0] + 1
        logger.info(
            f"Functional-group H detected on {atom['symbol']} "
            f"(CDXML id={atom['cdxml_id']}): 1-based Molecule index "
            f"{proton_index} (implicit_h_color={atom['implicit_h_color']})."
        )
        return proton_index

    def get_colored_proton_index(self, color_code=None, mode="acid"):
        """Identify the 1-based atom index of a site marked by colour.

        Supports two kinds of coloured protons:

        1. **Explicit H nodes** – a standalone ``<n Element="1">`` atom
           drawn in a different colour (e.g. a separate –H fragment).
        2. **Functional-group H** – a heavy-atom node whose label contains
           an "H" character in a *different* colour from the atom symbol
           (e.g. ``<s color="0">O</s><s color="4">H</s>`` for phenol –OH).
           The hydrogen is implicit (``NumHydrogens >= 1``) and only
           becomes an explicit atom after RDKit's ``Chem.AddHs``.

        Two detection modes are available:

        **Default mode** (``color_code=None``):
            Auto-detect: find the atom (explicit H *or* functional-group H)
            whose colour differs from the majority.

        **User-specified mode** (``color_code=<int>``):
            Find atoms / functional-group H spans matching that colour.

        Args:
            color_code (int | None): Colour-table index to look for.
                ``None`` triggers automatic detection.
            mode: ``"acid"`` (default) returns the uniquely coloured
                hydrogen. ``"base"`` returns the uniquely coloured
                heavy atom to protonate.

        Returns:
            int: 1-based index of the marked atom in the Molecule
                produced by ``CDXFile.get_molecules`` (i.e. after
                ``Chem.AddHs`` and 3-D embedding).

        Raises:
            ValueError: When the site cannot be identified unambiguously.
        """
        atoms = self.parse_cdxml_element_colors()

        if not atoms:
            raise ValueError(f"No atoms found in CDXML file: {self.filename}")

        if color_code is not None:
            return self._proton_index_by_color(atoms, color_code, mode=mode)
        return self._proton_index_auto(atoms, mode=mode)

    # ----- private helpers for get_colored_proton_index -----

    def _proton_index_by_color(self, atoms, color_code, mode="acid"):
        """Return the 1-based site index for a user-specified colour."""
        if mode == "base":
            heavy = [
                (i, a)
                for i, a in enumerate(atoms)
                if a["color"] == color_code and a["symbol"] != "H"
            ]
            if len(heavy) == 1:
                idx, atom = heavy[0]
                site = idx + 1
                logger.info(
                    f"User-specified colour: basic atom at 1-based index "
                    f"{site} (CDXML id={atom['cdxml_id']}, "
                    f"symbol={atom['symbol']}, color={color_code})."
                )
                return site
            if len(heavy) > 1:
                raise ValueError(
                    f"Multiple heavy atoms ({len(heavy)}) have "
                    f"color code {color_code}."
                )
            hydrogenish = [
                a
                for a in atoms
                if a["color"] == color_code and a["symbol"] == "H"
            ] or [a for a in atoms if a["implicit_h_color"] == color_code]
            if hydrogenish:
                raise ValueError(_PKB_COLORED_HYDROGEN_ERROR)
            raise ValueError(
                f"No atoms with color code {color_code} found in "
                f"{self.filename}."
            )

        # Case 1: explicit H node with that colour
        explicit_h = [
            (i, a)
            for i, a in enumerate(atoms)
            if a["color"] == color_code and a["symbol"] == "H"
        ]
        if len(explicit_h) == 1:
            idx, atom = explicit_h[0]
            proton_index = idx + 1
            logger.info(
                f"User-specified colour: explicit H at 1-based index "
                f"{proton_index} (CDXML id={atom['cdxml_id']}, "
                f"color={color_code})."
            )
            return proton_index
        if len(explicit_h) > 1:
            raise ValueError(
                f"Multiple explicit H atoms ({len(explicit_h)}) have "
                f"color code {color_code}."
            )

        # Case 2: functional-group implicit H with that colour
        fg_matches = [
            (i, a)
            for i, a in enumerate(atoms)
            if a["implicit_h_color"] == color_code
        ]
        if len(fg_matches) == 0:
            # Also check if any non-H atoms have that node colour (give
            # a clear message distinguishing "no match" from "not H")
            any_match = [a for a in atoms if a["color"] == color_code]
            if any_match:
                raise ValueError(
                    f"Atom(s) with color code {color_code} found, but "
                    f"none are hydrogen. Symbols: "
                    f"{[a['symbol'] for a in any_match]}."
                )
            raise ValueError(
                f"No atoms with color code {color_code} found in "
                f"{self.filename}."
            )
        if len(fg_matches) > 1:
            raise ValueError(
                f"Multiple functional-group H spans ({len(fg_matches)}) "
                f"have color code {color_code}."
            )

        cdxml_idx, atom = fg_matches[0]
        return self._resolve_implicit_h_index(atoms, cdxml_idx, atom)

    def _proton_index_auto(self, atoms, mode="acid"):
        """Auto-detect the uniquely coloured acid or base site."""
        # Gather *all* colours that appear – both node colours and
        # implicit-H span colours.
        color_counts = _cdxml_atom_color_counts(atoms)
        kind = (
            "basic atom to protonate" if mode == "base" else "proton to remove"
        )
        if len(color_counts) < 2:
            raise ValueError(
                "All atoms in the CDXML file share the same colour. "
                f"Cannot auto-detect the {kind}. "
                "Use -cl/--color-code to specify the colour explicitly, "
                "or use -pi/--proton-index."
            )

        majority_color = color_counts.most_common(1)[0][0]

        # Candidates: explicit H with a non-majority colour
        explicit_h = [
            (i, a)
            for i, a in enumerate(atoms)
            if a["symbol"] == "H" and a["color"] != majority_color
        ]

        # Candidates: functional-group implicit H whose span colour
        # differs from the majority
        fg_h = [
            (i, a)
            for i, a in enumerate(atoms)
            if a["implicit_h_color"] is not None
            and a["implicit_h_color"] != majority_color
        ]
        unique_heavy = [
            (i, a)
            for i, a in enumerate(atoms)
            if a["symbol"] != "H" and a["color"] != majority_color
        ]

        if mode == "base":
            if len(unique_heavy) == 1 and not explicit_h and not fg_h:
                idx, atom = unique_heavy[0]
                site = idx + 1
                logger.info(
                    f"Auto-detected basic atom at 1-based index "
                    f"{site} (CDXML id={atom['cdxml_id']}, "
                    f"symbol={atom['symbol']}, color={atom['color']})."
                )
                return site
            if len(unique_heavy) > 1:
                raise ValueError(
                    f"Multiple uniquely coloured heavy atoms found "
                    f"({len(unique_heavy)}). Cannot determine which atom "
                    f"to protonate. Use -cc/--color-code or "
                    f"-pi/--proton-index to specify explicitly."
                )
            if explicit_h or fg_h:
                raise ValueError(_PKB_COLORED_HYDROGEN_ERROR)
            raise ValueError(
                "No uniquely coloured atom found in the CDXML file."
            )

        total_candidates = len(explicit_h) + len(fg_h)

        if total_candidates == 0:
            # Collect all uniquely coloured atoms for the error message
            unique_atoms = [a for a in atoms if a["color"] != majority_color]
            if unique_atoms:
                symbols = [a["symbol"] for a in unique_atoms]
                raise ValueError(
                    f"Uniquely coloured atom(s) found but none are "
                    f"hydrogen (found: {symbols}). Only hydrogen atoms "
                    f"can be removed for pKa calculations."
                )
            raise ValueError(
                "No uniquely coloured atom found in the CDXML file."
            )

        if total_candidates > 1:
            raise ValueError(
                f"Multiple uniquely coloured hydrogen atoms found "
                f"({total_candidates}). Cannot determine which proton "
                f"to remove. Use -cl/--color-code or -pi/--proton-index "
                f"to specify explicitly."
            )

        # Exactly one candidate
        if explicit_h:
            idx, atom = explicit_h[0]
            proton_index = idx + 1
            logger.info(
                f"Auto-detected explicit H at 1-based index "
                f"{proton_index} (CDXML id={atom['cdxml_id']}, "
                f"color={atom['color']})."
            )
            return proton_index

        cdxml_idx, atom = fg_h[0]
        return self._resolve_implicit_h_index(atoms, cdxml_idx, atom)

    def _resolve_implicit_h_index(self, atoms, cdxml_idx, atom):
        """Map a functional-group implicit H to a 1-based Molecule index."""
        from rdkit import Chem

        num_h = atom.get("num_hydrogens")
        nested_ids = self._nested_deprotonatable_cdxml_ids(atom["cdxml_id"])
        if (num_h is None or num_h < 1) and not nested_ids:
            raise ValueError(
                f"Atom at CDXML position {cdxml_idx + 1} "
                f"(id={atom['cdxml_id']}, {atom['symbol']}) has a "
                f"coloured H in its label but NumHydrogens is "
                f"{num_h!r}. Cannot identify an implicit hydrogen."
            )

        rdkit_mols = list(
            Chem.MolsFromCDXMLFile(self.filename, removeHs=False)
        )
        if not rdkit_mols:
            raise ValueError(
                f"RDKit could not read {self.filename} for implicit-H "
                f"resolution."
            )

        for rdkit_mol in rdkit_mols:
            if rdkit_mol is None:
                continue
            rdkit_mol_h = Chem.AddHs(rdkit_mol)
            heavy_idx = self._rdkit_heavy_idx_for_implicit_h(
                rdkit_mol_h, atom, local_idx=cdxml_idx
            )
            if heavy_idx is None:
                continue
            try:
                return self._proton_index_from_rdkit_heavy(
                    rdkit_mol_h, heavy_idx, atom
                )
            except ValueError:
                continue

        raise ValueError(
            f"No hydrogen bonded to {atom['symbol']} "
            f"(CDXML id={atom['cdxml_id']}) after AddHs."
        )

    def get_colored_basic_atom_index(self, color_code=None):
        """Identify the 1-based index of a heavy atom marked by colour.

        Thin wrapper around :meth:`get_colored_proton_index` with
        ``mode="base"``.
        """
        return self.get_colored_proton_index(
            color_code=color_code, mode="base"
        )

    # ------------------------------------------------------------------
    # PKaMolecule factory
    # ------------------------------------------------------------------

    def get_pka_molecules(
        self,
        proton_index=None,
        color_code=None,
        index="-1",
        return_list=False,
    ):
        """Return :class:`PKaMolecule` instance(s) with the acidic proton resolved.

        Mirrors the signature and behaviour of :meth:`CDXFile.get_molecules`:

        * ``index='-1'`` (default) – returns a **single** ``PKaMolecule``.
        * ``index=':'``             – returns a **list** of ``PKaMolecule``.
        * ``return_list=True``      – always returns a list regardless of
          *index*.

        Proton detection priority:

        1. **User-supplied atom index** – ``proton_index`` (1-based).
        2. **User-supplied colour code** – ``color_code``; the matching
           coloured hydrogen is looked up via
           :meth:`get_colored_proton_index`.
        3. **Auto-detection** – :meth:`get_pka_molecules_auto` performs
           per-fragment colour analysis to identify the proton
           independently in each fragment.

        Args:
            proton_index (int | None): Explicit 1-based proton index.
                Skips all colour-based detection when provided.
            color_code (int | None): CDXML colour-table index.  Used
                only when *proton_index* is ``None``.
            index (str | slice): Fragment selector using the same
                1-based convention as :meth:`CDXFile.get_molecules`
                (``"-1"`` = last fragment, ``":"`` = all, ``"1"`` =
                first, etc.).  Defaults to ``"-1"``.
            return_list (bool): When ``True`` the return value is always
                a list even if *index* selects a single fragment.

        Returns:
            PKaMolecule | list[PKaMolecule]: A single ``PKaMolecule``
                when *index* selects one fragment and *return_list* is
                ``False``; a list otherwise.

        Raises:
            ValueError: If the proton cannot be identified or validated.
        """
        from chemsmart.io.molecules.structure import PKaMolecule

        # Build the full list of PKaMolecules for all fragments first.
        if proton_index is None and color_code is None:
            all_pka_mols = self.get_pka_molecules_auto()
        else:
            if proton_index is None:
                proton_index = self.get_colored_proton_index(
                    color_code=color_code,
                )
            all_pka_mols = [
                PKaMolecule(molecule=mol, proton_index=proton_index)
                for mol in self.molecules
            ]

        # Apply the same index semantics as CDXFile.get_molecules.
        if index == ":":
            selection = all_pka_mols
        else:
            parsed_index = (
                index
                if isinstance(index, slice)
                else string2index_1based(str(index))
            )
            selection = all_pka_mols[parsed_index]

        if return_list:
            return selection if isinstance(selection, list) else [selection]
        return selection

    # ------------------------------------------------------------------
    # Per-fragment colour parsing and auto-detection
    # ------------------------------------------------------------------

    def parse_cdxml_fragment_colors(self):
        """Parse per-fragment atom colour information from a CDXML file.

        Unlike :meth:`parse_cdxml_element_colors`, which returns a flat
        list across all fragments, this method returns a **list of lists**
        where each inner list contains the atom-colour dicts for one
        top-level ``<fragment>`` element.

        The per-fragment grouping preserves the 1-to-1 correspondence
        between fragments and the ``Molecule`` objects returned by
        :meth:`CDXFile.molecules` / ``Chem.MolsFromCDXMLFile``.

        Returns:
            list[list[dict]]: One sub-list per top-level fragment.
                Each dict has the same keys as
                :meth:`parse_cdxml_element_colors`.

        Raises:
            ValueError: If the file cannot be parsed as valid CDXML.
        """
        from chemsmart.utils.periodictable import PeriodicTable

        pt = PeriodicTable()
        root = self._parse_cdxml_root()

        fragments_atoms = []
        for fragment in self._iter_top_level_fragments(root):
            fragment_atoms = self._parse_fragment_nodes(fragment, pt)
            if fragment_atoms:
                fragments_atoms.append(fragment_atoms)

        return fragments_atoms

    def _detect_proton_in_fragment(
        self, atoms, fragment_index=None, mode="acid"
    ):
        """Auto-detect the uniquely coloured site within a single fragment.

        The logic mirrors :meth:`_proton_index_auto` but operates on a
        fragment-local atom list and returns a **fragment-local** 0-based
        index (the offset within the fragment's atom list before
        ``AddHs``).

        For implicit / functional-group hydrogens the method returns a
        tuple ``(heavy_atom_local_idx, atom_dict)`` so the caller can
        resolve the final Molecule index via RDKit ``AddHs``.

        Args:
            atoms: List of atom dicts for one fragment (from
                :meth:`parse_cdxml_fragment_colors`).
            fragment_index: Optional fragment number for error messages.
            mode: ``"acid"`` (coloured hydrogen) or ``"base"`` (coloured
                heavy atom).

        Returns:
            dict: With keys:
                * ``type`` – ``"explicit"``, ``"implicit"``, or ``"heavy"``
                * ``local_idx`` – 0-based position in the fragment's
                  atom list (before ``AddHs``)
                * ``atom`` – the atom dict

        Raises:
            ValueError: When the site cannot be identified.
        """
        frag_info = (
            f" (fragment {fragment_index})"
            if fragment_index is not None
            else ""
        )

        color_counts = _cdxml_atom_color_counts(atoms)
        kind = (
            "basic atom to protonate" if mode == "base" else "proton to remove"
        )
        if len(color_counts) < 2:
            raise ValueError(
                f"All atoms in fragment{frag_info} share the same colour. "
                f"Cannot auto-detect the {kind}."
            )

        majority_color = color_counts.most_common(1)[0][0]

        explicit_h = [
            (i, a)
            for i, a in enumerate(atoms)
            if a["symbol"] == "H" and a["color"] != majority_color
        ]
        # Functional-group implicit H on *heavy* atoms only — skip
        # explicit H nodes to avoid double-counting.
        fg_h = [
            (i, a)
            for i, a in enumerate(atoms)
            if a["symbol"] != "H"
            and a["implicit_h_color"] is not None
            and a["implicit_h_color"] != majority_color
        ]
        unique_heavy = [
            (i, a)
            for i, a in enumerate(atoms)
            if a["symbol"] != "H" and a["color"] != majority_color
        ]

        if mode == "base":
            if len(unique_heavy) == 1 and not explicit_h and not fg_h:
                idx, atom = unique_heavy[0]
                return {"type": "heavy", "local_idx": idx, "atom": atom}
            if len(unique_heavy) > 1:
                raise ValueError(
                    f"Multiple uniquely coloured heavy atoms "
                    f"({len(unique_heavy)}) in fragment{frag_info}. "
                    f"Cannot determine which atom to protonate."
                )
            if explicit_h or fg_h:
                raise ValueError(_PKB_COLORED_HYDROGEN_ERROR)
            raise ValueError(
                f"No uniquely coloured atom found in fragment{frag_info}."
            )

        total = len(explicit_h) + len(fg_h)

        if total == 0:
            unique = [a for a in atoms if a["color"] != majority_color]
            if unique:
                raise ValueError(
                    f"Uniquely coloured atom(s) found in fragment{frag_info} "
                    f"but none are hydrogen (found: "
                    f"{[a['symbol'] for a in unique]})."
                )
            raise ValueError(
                f"No uniquely coloured atom found in fragment{frag_info}."
            )

        if total > 1:
            raise ValueError(
                f"Multiple uniquely coloured H atoms ({total}) in "
                f"fragment{frag_info}. Cannot determine which proton "
                f"to remove."
            )

        if explicit_h:
            idx, atom = explicit_h[0]
            return {"type": "explicit", "local_idx": idx, "atom": atom}

        idx, atom = fg_h[0]
        return {"type": "implicit", "local_idx": idx, "atom": atom}

    def _detect_basic_atom_in_fragment(self, atoms, fragment_index=None):
        """Thin wrapper around :meth:`_detect_proton_in_fragment` for ``--pkb``."""
        return self._detect_proton_in_fragment(
            atoms, fragment_index=fragment_index, mode="base"
        )

    def _rdkit_mols_with_hs(self):
        """Return RDKit molecules from this file after ``Chem.AddHs``."""
        from rdkit import Chem

        rdkit_mols = list(
            Chem.MolsFromCDXMLFile(self.filename, removeHs=False)
        )
        result = []
        for rdkit_mol in rdkit_mols:
            if rdkit_mol is not None:
                result.append(Chem.AddHs(rdkit_mol))
            else:
                result.append(None)
        return result

    def _site_index_from_fragment_detection(
        self, detection, rdkit_mol_h, fragment_index
    ):
        """Map a fragment colour detection to a 1-based Molecule index."""
        if rdkit_mol_h is None:
            raise ValueError(
                f"RDKit molecule for fragment {fragment_index} is "
                f"None; cannot resolve site index."
            )
        if detection["type"] != "implicit":
            heavy_idx = self._rdkit_atom_idx_by_cdxml_id(
                rdkit_mol_h, detection["atom"]["cdxml_id"]
            )
            if heavy_idx is None:
                heavy_idx = detection["local_idx"]
            return heavy_idx + 1
        heavy_idx = self._rdkit_heavy_idx_for_implicit_h(
            rdkit_mol_h,
            detection["atom"],
            local_idx=detection["local_idx"],
        )
        if heavy_idx is None:
            raise ValueError(
                f"No hydrogen bonded to "
                f"{detection['atom']['symbol']} "
                f"(CDXML id={detection['atom']['cdxml_id']}) "
                f"in fragment {fragment_index} after AddHs."
            )
        return self._proton_index_from_rdkit_heavy(
            rdkit_mol_h, heavy_idx, detection["atom"]
        )

    def _iter_colored_sites(self, color_code=None, mode="acid"):
        """Yield ``(molecule, 1-based site index)`` for each CDXML fragment."""
        molecules = list(self.molecules)
        if color_code is not None:
            site = self.get_colored_proton_index(
                color_code=color_code, mode=mode
            )
            for mol in molecules:
                yield mol, site
            return

        fragments_atoms = self.parse_cdxml_fragment_colors()
        if len(fragments_atoms) != len(molecules):
            logger.warning(
                f"Fragment count ({len(fragments_atoms)}) differs from "
                f"molecule count ({len(molecules)}). Falling back to "
                f"global colour detection."
            )
            site = self.get_colored_proton_index(mode=mode)
            for mol in molecules:
                yield mol, site
            return

        rdkit_mols_h = self._rdkit_mols_with_hs()
        for frag_idx, (frag_atoms, mol, rdkit_mol_h) in enumerate(
            zip(fragments_atoms, molecules, rdkit_mols_h)
        ):
            detection = self._detect_proton_in_fragment(
                frag_atoms, fragment_index=frag_idx + 1, mode=mode
            )
            site = self._site_index_from_fragment_detection(
                detection, rdkit_mol_h, frag_idx + 1
            )
            logger.info(
                f"Fragment {frag_idx + 1}: detected site at 1-based "
                f"index {site} "
                f"(type={detection['type']}, "
                f"cdxml_id={detection['atom']['cdxml_id']})."
            )
            yield mol, site

    def get_pka_molecules_auto(self, mode="acid"):
        """Per-fragment colour auto-detection → list of PKaMolecule.

        For each top-level ``<fragment>`` in the CDXML file:

        1. Parse atom colours within that fragment.
        2. Identify the uniquely coloured site (acidic H, or with
           ``mode="base"`` the heavy atom to protonate).
        3. Map the site to the 1-based index in the ``Molecule``
           produced after ``Chem.AddHs`` + embedding.
        4. Wrap as a :class:`PKaMolecule`. For ``mode="base"`` the
           molecule is protonated first so ``proton_index`` is the
           added hydrogen.

        Falls back to :meth:`get_colored_proton_index` if per-fragment
        parsing produces a fragment/molecule count mismatch.

        Returns:
            list[PKaMolecule]: One ``PKaMolecule`` per fragment.

        Raises:
            ValueError: If detection fails for any fragment.
        """
        from chemsmart.io.molecules.structure import PKaMolecule

        if mode == "base":
            return [
                PKaMolecule.add_proton_at_atom(mol, site)
                for mol, site in self._iter_colored_sites(mode="base")
            ]
        return [
            PKaMolecule(molecule=mol, proton_index=site)
            for mol, site in self._iter_colored_sites(mode="acid")
        ]

    def get_pkb_molecules_auto(self):
        """Per-fragment basic-atom auto-detection → protonated PKaMolecules."""
        return self.get_pka_molecules_auto(mode="base")
