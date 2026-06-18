# ORCA Ionic-Crystal QM/MM `--geninput`

> **User documentation:** see [docs/orca-ionic-crystal-qmmm-geninput.md](docs/orca-ionic-crystal-qmmm-geninput.md)
> for the full workflow, architecture, and NaCl example command.

This file is retained as a short pointer. The implemented design uses a
**Python prejob hook** inside `chemsmart_run_{label}.py`, not CrystalPrep
commands in `chemsmart_sub_{label}.sh`.
