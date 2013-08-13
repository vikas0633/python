
import threading
import os

semaphore = threading.Semaphore(4)

def run_command(cmd):
    with semaphore:
        os.system(cmd)

for i in range(8):
    threading.Thread(target=run_command, args=("sleep 10", )).start()
