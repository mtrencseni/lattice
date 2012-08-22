cd latex
for /f %%f in ('dir /b ..\scripts\*.gnu') do wgnuplot ..\scripts\%%f
cd ..