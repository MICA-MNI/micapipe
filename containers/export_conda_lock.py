#!/usr/bin/env python3
"""Export an explicit URL and checksum lock from installed conda metadata."""
import json
import sys
from pathlib import Path
prefix=Path(sys.argv[1])
records=[json.loads(p.read_text()) for p in sorted((prefix/'conda-meta').glob('*.json'))]
if not records:
    raise SystemExit('No conda package records found')
lines=['# platform: linux-64','@EXPLICIT']
for rec in records:
    url=rec['url']
    checksum=rec.get('md5')
    if checksum:
        url += '#'+checksum
    lines.append(url)
Path(sys.argv[2]).write_text('\n'.join(lines)+'\n')
print(f'{len(records)} packages locked')
