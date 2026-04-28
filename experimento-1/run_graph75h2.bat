@echo off
setlocal enabledelayedexpansion
set EXE=.\sgap-heuristic.exe
set CONFIG=.\greedy_only.ini
set SAIDA=.\resultados_graph75h2.csv

for %%G in (.\graph-75\*.gfa) do (
    for %%N in (5 10 15 20 25 30 35 40 45 50 55 60 65 70 75) do (
        for %%S in (.\seq-%%N\*.fastq) do (
            echo %%~nxG ^<- %%~nxS
            "%EXE%" "%%G" "%%S" "%CONFIG%" "%SAIDA%"
            if errorlevel 1 (
                echo ERRO na execucao
                exit /b 1
            )
        )
    )
)
echo Concluido. Resultados em %SAIDA%
endlocal
