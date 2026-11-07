"""Regenerate every paper figure from the data:  python paper/scripts/make_all.py

Fig. 1 re-solves two modes and needs the femwell environment (FEMWELL_PYTHON, default the project venv); the others
run in the main environment.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FEMWELL = os.environ.get('FEMWELL_PYTHON', 'C:/Users/unmol/.venvs/femwell/Scripts/python.exe')
SCRIPTS = [('fig1_sections.py', FEMWELL), ('fig2_ratio.py', sys.executable), ('fig3_stages.py', sys.executable),
           ('fig4_montecarlo.py', sys.executable), ('fig5_slab.py', sys.executable),
           ('fig6_variation.py', sys.executable), ('fig7_noise.py', sys.executable)]

if __name__ == '__main__':
    env = dict(os.environ, PYTHONIOENCODING='utf-8')
    for script, py in SCRIPTS:
        print(f'--> {script}', flush=True)
        subprocess.run([py, os.path.join(HERE, script)], cwd=HERE, env=env, check=True)
