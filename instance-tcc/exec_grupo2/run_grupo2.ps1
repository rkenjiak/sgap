# Executa beam, brkga, jain e pli sobre as instancias do grupo 2 (pli so no no raiz).
#
# Uso (a partir de qualquer pasta):
#   powershell -ExecutionPolicy Bypass -File run_grupo2.ps1
#   powershell -ExecutionPolicy Bypass -File run_grupo2.ps1 -Algoritmos beam,jain
#   powershell -ExecutionPolicy Bypass -File run_grupo2.ps1 -Algoritmos brkga -SementesBrkga 1..3
#
# Pareamento: para cada pasta de tamanho (ex.: 125s_250v), a sequencia
# <N>s_a4_inst<i>.fastq eh alinhada aos 5 grafos <M>v_a4_d*_inst<i>.gfa.
#
# BRKGA: eh deterministico para uma semente, entao cada run usa uma
# semente diferente (padrao 1..10). Para cada semente o script gera
# config/brkga_sementes/brkga_s<semente>.ini a partir de config/brkga.ini
# (trocando so brkga_seed). Todas as runs vao para resultados/brkga.csv;
# a semente aparece na coluna variant (ex.: vp_p256_g500_s3).
#
# Saida: resultados/<algoritmo>.csv (uma linha por grafo x sequencia,
# e por semente no BRKGA). Execucoes ja presentes no CSV sao puladas,
# entao o script pode ser interrompido e retomado.

param(
  [string[]]$Algoritmos    = @("beam", "brkga", "jain", "pli"),
  [string[]]$SementesBrkga = @("1..10")
)

$ErrorActionPreference = "Stop"

# Com -File, "beam,jain" e "1..10" chegam como strings
$Algoritmos = $Algoritmos | ForEach-Object { $_ -split "," } | Where-Object { $_ } | ForEach-Object { $_.Trim().ToLower() }
$SementesBrkga = $SementesBrkga | ForEach-Object { $_ -split "," } | Where-Object { $_ } | ForEach-Object {
  $s = $_.Trim()
  if ($s -match '^(\d+)\.\.(\d+)$') { [int]$Matches[1]..[int]$Matches[2] } else { [int]$s }
}

$Grupo      = "grupo2"
$Raiz       = Split-Path -Parent $PSScriptRoot
$BinDir     = Join-Path $Raiz "sgap\x64\Release"
$InstDir    = Join-Path $Raiz $Grupo
$ConfigDir  = Join-Path $PSScriptRoot "config"
$ResultDir  = Join-Path $PSScriptRoot "resultados"

$Binario = @{
  beam  = Join-Path $BinDir "sgap-heuristic.exe"
  brkga = Join-Path $BinDir "sgap-heuristic.exe"
  jain  = Join-Path $BinDir "sgap-heuristic.exe"
  pli   = Join-Path $BinDir "sgap-solver.exe"
}

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

Write-Host "[$Grupo] $($Pares.Count) pares grafo x sequencia; algoritmos: $($Algoritmos -join ', ')"
if ($Algoritmos -contains "brkga") { Write-Host "[$Grupo] sementes do brkga: $($SementesBrkga -join ', ')" }

# Monta a lista de rodadas: uma por algoritmo, e uma por semente no BRKGA.
# Sufixo: final esperado da coluna variant (usado para saber o que ja rodou).
$Rodadas = @()
foreach ($alg in $Algoritmos) {
  if (-not $Binario.ContainsKey($alg)) { throw "algoritmo desconhecido: $alg" }
  $base = Join-Path $ConfigDir "$alg.ini"

  if ($alg -eq "brkga") {
    $dirSem = Join-Path $ConfigDir "brkga_sementes"
    New-Item -ItemType Directory -Force $dirSem | Out-Null
    $linhas = Get-Content $base | Where-Object { $_ -notmatch '^\s*brkga_seed\s*=' }
    foreach ($s in $SementesBrkga) {
      $cfgSem = Join-Path $dirSem "brkga_s$s.ini"
      Set-Content $cfgSem -Encoding ascii -Value (@("# GERADO por run_grupo2.ps1 a partir de ../brkga.ini - nao editar") + $linhas + "brkga_seed = $s")
      $Rodadas += [pscustomobject]@{ Alg = $alg; Nome = "brkga (semente $s)"; Cfg = $cfgSem; Sufixo = "_s$s" }
    }
  } else {
    $Rodadas += [pscustomobject]@{ Alg = $alg; Nome = $alg; Cfg = $base; Sufixo = "" }
  }
}

foreach ($r in $Rodadas) {
  $exe = $Binario[$r.Alg]
  $csv = Join-Path $ResultDir "$($r.Alg).csv"
  if (-not (Test-Path $exe)) { throw "executavel nao encontrado: $exe" }

  # Pares ja executados nesta rodada (graph,sequence + final do variant)
  $feitos = @{}
  if (Test-Path $csv) {
    Get-Content $csv | Select-Object -Skip 1 | ForEach-Object {
      $c = $_.Split(',')
      if ($c.Count -gt 6 -and $c[6].EndsWith($r.Sufixo)) { $feitos["$($c[0]),$($c[1])"] = $true }
    }
  }

  Write-Host "`n===== $($r.Nome) ($exe) =====" -ForegroundColor Cyan
  $k = 0
  $t0 = Get-Date
  foreach ($p in $Pares) {
    $k++
    $chave = "$($p.Gfa.BaseName),$($p.Fastq.BaseName)"
    if ($feitos.ContainsKey($chave)) { continue }

    Write-Host ("[{0}/{1}] {2} {3}" -f $k, $Pares.Count, $r.Nome, $chave)
    & $exe $p.Gfa.FullName $p.Fastq.FullName $r.Cfg $csv
    if ($LASTEXITCODE -ne 0) {
      Write-Warning "falha (exit $LASTEXITCODE) em $($r.Nome) $chave"
    }
  }
  Write-Host ("$($r.Nome) concluido em {0:N1} min" -f ((Get-Date) - $t0).TotalMinutes) -ForegroundColor Green
}
