# fix_venv.ps1 — このプロジェクトの .venv を「現マシンのローカル Python 3.13」に向け直す。
#
# 用途: PC を変えた直後など、.venv\Scripts\python.exe を叩くと
#   did not find executable at 'C:\Users\<別ユーザー>\...\python.exe'
# で起動失敗する時に 1 回実行する。
#
# 仕組み: venv はベース Python の絶対パスを pyvenv.cfg に焼き込む (= マシン固有)。
# .venv フォルダを PC 間で持ち回ると別ユーザーのパスになり壊れる。本スクリプトは
# py ランチャ (py -3.13) で現マシンの Python 3.13 を検出し pyvenv.cfg を書き換える。
# site-packages (optuna 等) は一切触らない・完全に可逆。cp313 ホイールは 3.13.x 共通
# なので、3.13 系である限りパッチ違い (3.13.0 vs 3.13.13 等) は ABI 互換で問題なし。
#
# 使い方:  PS> .\fix_venv.ps1

$ErrorActionPreference = 'Stop'

$venvCfg = Join-Path $PSScriptRoot '.venv\pyvenv.cfg'
$venvPy  = Join-Path $PSScriptRoot '.venv\Scripts\python.exe'
if (-not (Test-Path $venvCfg)) {
    Write-Host "[fix_venv] .venv\pyvenv.cfg が見つかりません ($venvCfg)。" -ForegroundColor Red
    Write-Host "[fix_venv] .venv 未作成かも。 py -3.13 -m venv .venv で作成してください。" -ForegroundColor Red
    exit 1
}

# --- 現マシンの Python 3.13 を検出 (パス/ユーザー名に依存しない) ---
# .venv の site-packages は cp313 ホイールのため 3.13 系であることが必須。
# py -3.13 が無ければ 3.13 未導入 → 明示エラー (3.12/3.14 を勝手に使わない)。
$pyExe = ''
try {
    $pyExe = ((& py -3.13 -c "import sys; print(sys.executable)") -join '').Trim()
} catch {
    $pyExe = ''
}
if ([string]::IsNullOrWhiteSpace($pyExe) -or -not (Test-Path $pyExe)) {
    Write-Host "[fix_venv] Python 3.13 が見つかりません (py -3.13 が失敗)。" -ForegroundColor Red
    Write-Host "[fix_venv] このマシンに Python 3.13.x を入れてから再実行してください。" -ForegroundColor Red
    exit 1
}
$pyVer   = ((& $pyExe -c "import sys; print('.'.join(map(str, sys.version_info[:3])))") -join '').Trim()
$pyHome  = Split-Path $pyExe -Parent
$venvDir = Join-Path $PSScriptRoot '.venv'

Write-Host "[fix_venv] 検出: Python $pyVer  ($pyExe)" -ForegroundColor Cyan

# --- pyvenv.cfg を書き換え (include-system-site-packages は false 維持) ---
$content = @"
home = $pyHome
include-system-site-packages = false
version = $pyVer
executable = $pyExe
command = $pyExe -m venv $venvDir
"@
# BOM 無し UTF-8 で書く (venv ランチャが先頭 BOM で home 行を読み損なうのを防ぐ)
[System.IO.File]::WriteAllText($venvCfg, $content, (New-Object System.Text.UTF8Encoding($false)))
Write-Host "[fix_venv] pyvenv.cfg を更新しました (home=$pyHome)" -ForegroundColor Green

# --- 検証: venv python が起動し主要依存が import できるか ---
try {
    $check = ((& $venvPy -c "import sys, optuna; print('python ' + sys.version.split()[0] + ' / optuna ' + optuna.__version__)") -join '').Trim()
    Write-Host "[fix_venv] 検証 OK: $check" -ForegroundColor Green
    Write-Host "[fix_venv] 完了。 .\.venv\Scripts\python.exe main.py で実行できます。" -ForegroundColor Green
} catch {
    Write-Host "[fix_venv] 検証で venv python の起動/import に失敗しました:" -ForegroundColor Red
    Write-Host $_.Exception.Message -ForegroundColor Red
    Write-Host "[fix_venv] site-packages 破損の可能性 → py -3.13 -m venv --clear .venv で作り直し後、依存を再インストールしてください。" -ForegroundColor Red
    exit 1
}
