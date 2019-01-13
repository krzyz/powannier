"set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)" > newFile.txt
Get-Content CMakeLists.txt >> newFile.txt
Get-Content newFile.txt > CMakeLists.txt
(Get-Content -Encoding UTF8 CMakeLists.txt) -replace 'cubature m', 'cubature' | Out-File -Encoding UTF8 CMakeLists.txt
(Get-Content -Encoding UTF8 CMakeLists.txt) -replace 'cubature SHARED', 'cubature STATIC' | Out-File -Encoding UTF8 CMakeLists.txt
