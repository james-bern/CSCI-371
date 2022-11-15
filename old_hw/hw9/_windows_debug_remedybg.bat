@echo off
tasklist | find /i "remedybg.exe" && start /b remedybg.exe start-debugging || start /b remedybg.exe executable.exe && timeout /t 1 /nobreak && start /b remedybg.exe start-debugging
