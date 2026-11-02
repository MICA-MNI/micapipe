#!/usr/bin/env python3
"""Download model assets from a pinned LAMAR-Models revision."""
import ast,hashlib,json,shutil,urllib.request
from pathlib import Path
revision='5f43003b01e56d2754b4c0e23a98f38d1d134de3'
root=Path(__file__).resolve().parents[1]
workspace=root.parent
handler=workspace/'downloads/vendor-src/lamareg/lamareg/SynthSeg/model_handler.py'
tree=ast.parse(handler.read_text())
urls=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='MODEL_URLS' for t in n.targets))
destination=workspace/'downloads/lamareg-models'
destination.mkdir(exist_ok=True)
checksums={}
for name,sources in urls.items():
    target=destination/name
    if not target.exists():
        if isinstance(sources,str):sources=[sources]
        with target.with_suffix('.partial').open('wb') as output:
            for url in sources:
                pinned=url.replace('/refs/heads/main/',f'/{revision}/')
                print(pinned,flush=True)
                with urllib.request.urlopen(pinned,timeout=120) as stream:
                    shutil.copyfileobj(stream,output,16*1024*1024)
        target.with_suffix('.partial').replace(target)
    h=hashlib.sha256()
    with target.open('rb') as stream:
        for block in iter(lambda:stream.read(16*1024*1024),b''):h.update(block)
    checksums[name]=h.hexdigest()
(root/'containers/lamareg-models.json').write_text(json.dumps({'revision':revision,'sha256':checksums},indent=2)+'\n')
