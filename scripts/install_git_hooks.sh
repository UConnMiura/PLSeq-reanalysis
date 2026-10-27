#!/usr/bin/env bash
set -euo pipefail

root=$(git rev-parse --show-toplevel)
hook="$root/.git/hooks/pre-commit"

cat > "$hook" <<'HOOK'
#!/usr/bin/env bash
set -euo pipefail
root=$(git rev-parse --show-toplevel)
python3 "$root/tests/check_repository_files.py" --staged
HOOK

chmod +x "$hook"
printf 'Installed repository safety hook: %s\n' "$hook"

