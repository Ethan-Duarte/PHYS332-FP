# V1.5.1: macOS / Windows cross-platform package

This is a platform-support update to V1.5. It does **not** change force algorithms, galaxy initial conditions, numerical defaults, the 3D viewer, or saved-trajectory formats.

New:
- `start_gui.command` — Finder double-click entry point for macOS.
- `start_gui.sh` — creates `.venv/bin/python`, installs ordinary dependencies, launches same local browser GUI, keeps errors visible.
- `scripts/install_agama_macos.sh` — optional source installation with Apple compiler and Homebrew headers/libraries.
- `docs/MACOS_SETUP.md` — installation, Intel/Apple Silicon guidance, troubleshooting, and moving runs between OSes.
- Platform-aware AGAMA status text in the GUI, expanded documentation, and shell packaging tests.

Verification: 36 Python tests passed in a Linux environment; local GUI answered HTTP 200. The shell scripts pass `bash -n`; executable bits are included in the release ZIP. The Windows launcher is unchanged. **No native macOS GUI session or compiled macOS AGAMA library was available for testing.** AGAMA integration remains a candidate until its real backend and isolated equilibrium validation succeed.
