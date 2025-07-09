# Requires typst installed.
cd "$(dirname "${BASH_SOURCE[0]}")"
typst compile binder.typ ../print_format/binder.pdf
