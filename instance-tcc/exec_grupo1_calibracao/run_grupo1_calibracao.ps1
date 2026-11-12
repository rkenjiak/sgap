# Calibracao: executa apenas o pli sobre as instancias do grupo 1 com a mesma
# configuracao do grupo 2 (pli_node_limit = 0, so o no raiz). Comparando com
# exec_grupo1/resultados/pli.csv (otimo, sem limite de nos), da para medir o
# quanto o no raiz se afasta do otimo antes de olhar o grupo 2.
#
# Uso (a partir de qualquer pasta):
#   powershell -ExecutionPolicy Bypass -File run_grupo1_calibracao.ps1
#
# Pareamento: para cada pasta de tamanho (ex.: 25s_50v), a sequencia
# <N>s_a4_inst<i>.fastq eh alinhada aos 5 grafos <M>v_a4_d*_inst<i>.gfa.
#
# Saida: resultados/pli.csv (uma linha por grafo x sequencia). Execucoes ja
# presentes no CSV sao puladas, entao o script pode ser interrompido e retomado.

$ErrorActionPreference = "Stop"

$Grupo      = "grupo1"
$Raiz       = Split-Path -Parent $PSScriptRoot
$Exe        = Join-Path $Raiz "sgap\x64\Release\sgap-solver.exe"
$InstDir    = Join-Path $Raiz $Grupo
$Cfg        = Join-Path $PSScriptRoot "config\pli.ini"
$ResultDir  = Join-Path $PSScriptRoot "resultados"
$Csv        = Join-Path $ResultDir "pli.csv"

if (-not (Test-Path $Exe)) { throw "executavel nao encontrado: $Exe" }
New-Item -ItemType Directory -Force $ResultDir | Out-Null

# Lista de pares (gfa, fastq) do grupo
$Pares = @()
foreach ($dir in Get-ChildItem $InstDir -Directory | Sort-Object { [int]($_.Name -replace 's_.*', '') }) {
  foreach ($fq in Get-ChildItem $dir.FullName -Filter "*.fastq") {
    $inst = [regex]::Match($fq.BaseName, '_inst(\d+)$').Groups[1].Value
    $gfas = Get-ChildItem $dir.FullName -Filter "*_inst$inst.gfa" | Sort-Object Name
    foreach ($gfa in $gfas) {
      $Pares += [pscustomobject]@{ Gfa = $gfa; Fastq = $fq; Inst = [int]$inst; Dir = $dir.Name }
    }
  }
}
$Pares = $Pares | Sort-Object { [int]($_.Dir -replace 's_.*', '') }, Inst, { $_.Gfa.Name }

Write-Host "[$Grupo calibracao] $($Pares.Count) pares grafo x sequencia; pli no raiz"

# Pares ja executados (graph,sequence)
$feitos = @{}
if (Test-Path $Csv) {
  Get-Content $Csv | Select-Object -Skip 1 | ForEach-Object {
    $c = $_.Split(',')
    if ($c.Count -gt 1) { $feitos["$($c[0]),$($c[1])"] = $true }
  }
}

Write-Host "`n===== pli no raiz ($Exe) =====" -ForegroundColor Cyan
$k = 0
$t0 = Get-Date
foreach ($p in $Pares) {
  $k++
  $chave = "$($p.Gfa.BaseName),$($p.Fastq.BaseName)"
  if ($feitos.ContainsKey($chave)) { continue }

  Write-Host ("[{0}/{1}] pli {2}" -f $k, $Pares.Count, $chave)
  & $Exe $p.Gfa.FullName $p.Fastq.FullName $Cfg $Csv
  if ($LASTEXITCODE -ne 0) {
    Write-Warning "falha (exit $LASTEXITCODE) em pli $chave"
  }
}
Write-Host ("pli concluido em {0:N1} min" -f ((Get-Date) - $t0).TotalMinutes) -ForegroundColor Green
