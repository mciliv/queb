#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

detect_rc_file() {
  shell_path="${SHELL:-}"
  if [ -n "${ZDOTDIR:-}" ] && [ -f "${ZDOTDIR}/.zshrc" ]; then
    echo "${ZDOTDIR}/.zshrc"; return
  fi
  if [ "${shell_path##*/}" = "zsh" ]; then
    echo "${HOME}/.zshrc"; return
  fi
  if [ -f "${HOME}/.bashrc" ]; then
    echo "${HOME}/.bashrc"; return
  fi
  echo "${HOME}/.bash_profile"
}

rc_file="$(detect_rc_file)"
touch "$rc_file"
backup="${rc_file}.bak.$(date +%Y%m%d_%H%M%S)"
cp "$rc_file" "$backup"

tmp_file="$(mktemp)"
awk 'BEGIN{skip=0} /# mol aliases start/{skip=1; next} /# mol aliases end/{skip=0; next} skip==0{print}' "$rc_file" > "$tmp_file"
mv "$tmp_file" "$rc_file"

{
  echo
  echo "# mol aliases start"
  echo "MOL_ROOT=\"$project_root\""
  echo "mol_is_in_mol_dir() { case \"\$PWD/\" in \"$project_root/*\") return 0;; \*) return 1;; esac }"
  echo "mol_activate() { case :\"\$PATH\": in *:\"$project_root/scripts\":*) ;; *) export PATH=\"$project_root/scripts:\$PATH\" ;; esac; unalias dev 2>/dev/null || true; alias d=\"dev\" 2>/dev/null || true; hash -r 2>/dev/null || true; rehash 2>/dev/null || true; }"
  echo "mol_deactivate() { local p=:\"\$PATH\":; p=\"\${p//:$project_root\\/scripts:/:}\"; p=\"\${p#:}\"; p=\"\${p%:}\"; export PATH=\"\$p\"; unalias d 2>/dev/null || true; }"
  echo "mol_chpwd() { if mol_is_in_mol_dir; then mol_activate; else mol_deactivate; fi }"
  echo "# zsh hook"
  echo "if [ -n \"\$ZSH_VERSION\" ]; then autoload -U add-zsh-hook 2>/dev/null || true; add-zsh-hook chpwd mol_chpwd 2>/dev/null || true; fi"
  echo "# bash fallback"
  echo "if [ -n \"\$BASH_VERSION\" ]; then case \"\$PROMPT_COMMAND\" in *mol_chpwd*) ;; *) PROMPT_COMMAND=\"mol_chpwd;\$PROMPT_COMMAND\" ;; esac; fi"
  echo "# run once for current shell"
  echo "mol_chpwd"
  echo "# mol aliases end"
} >> "$rc_file"

echo "Aliases installed to $rc_file"
echo "Reload with: source \"$rc_file\""
