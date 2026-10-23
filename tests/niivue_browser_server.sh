#!/bin/bash
set -euo pipefail
cd "${1:?Pass the isolated server workspace}"
export UV_CACHE_DIR="$PWD/cache/uv" UV_PYTHON_INSTALL_DIR="$PWD/environments/python" UV_LINK_MODE=copy
test -x environments/browser/bin/python || tools/uv venv --python 3.11.14 environments/browser
tools/uv pip install --python environments/browser/bin/python playwright==1.55.0
export PLAYWRIGHT_BROWSERS_PATH="$PWD/cache/playwright"
environments/browser/bin/python -m playwright install chromium
environments/micapipe/bin/python - <<'PY'
import numpy as np,nibabel as nib
from pathlib import Path
out=Path('validation/niivue')
out.mkdir(parents=True,exist_ok=True)
data=np.zeros((32,32,32),dtype=np.float32)
data[8:24,8:24,8:24]=1
nib.save(nib.Nifti1Image(data,np.eye(4)),out/'phantom.nii.gz')
PY
environments/micapipe/bin/python source/functions/niivue_report.py --volume validation/niivue/phantom.nii.gz --output validation/niivue/viewer-v2
environments/browser/bin/python - <<'PY'
import functools,http.server,json,threading
from pathlib import Path
from playwright.sync_api import sync_playwright
handler=functools.partial(http.server.SimpleHTTPRequestHandler,directory='validation/niivue/viewer-v2')
server=http.server.ThreadingHTTPServer(('127.0.0.1',0),handler)
threading.Thread(target=server.serve_forever,daemon=True).start()
errors=[]; remote=[]
try:
 with sync_playwright() as p:
  browser=p.chromium.launch(headless=True,args=['--no-sandbox','--use-angle=swiftshader','--enable-unsafe-swiftshader'])
  page=browser.new_page(viewport={'width':1200,'height':800})
  page.on('pageerror',lambda err:errors.append(str(err)))
  def route(request):
   if request.request.url.startswith('http://127.0.0.1:'): request.continue_()
   else: remote.append(request.request.url);request.abort()
  page.route('**/*',route)
  page.goto(f'http://127.0.0.1:{server.server_port}/')
  page.wait_for_function("document.body.dataset.loaded === 'true'",timeout=60000)
  page.screenshot(path='evidence/niivue.png')
  Path('evidence/niivue-browser.json').write_text(json.dumps({'loaded':True,'page_errors':errors,'external_requests':remote},indent=2))
  assert not errors, errors
  assert not remote, remote
  browser.close()
finally: server.shutdown()
PY
