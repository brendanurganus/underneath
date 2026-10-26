"""Assemble site/index.html from src/page.css, src/body.html and src/app.js.

Usage: python3 assemble.py ROOT
"""
import os
import sys

root = sys.argv[1]
src = os.path.join(root, 'src')
css = open(os.path.join(src, 'page.css'), encoding='utf-8').read()
body = open(os.path.join(src, 'body.html'), encoding='utf-8').read()
js = open(os.path.join(src, 'app.js'), encoding='utf-8').read()
head = '''<title>Underneath</title>
<meta name="claude:viewport-fit-sides" content="cover">
<meta name="description" content="Every M4.5+ earthquake of the past year on a 3D globe, sliced open to the core, with real seismograms.">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Fredoka:wght@400;500;600&family=Nunito:wght@400;600;700;800&display=swap">
'''
importmap = '''<script type="importmap">
{ "imports": { "three": "https://cdn.jsdelivr.net/npm/three@0.170.0/build/three.module.js", "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.170.0/examples/jsm/" } }
</script>
'''
html = head + '<style>\n' + css + '\n</style>\n' + body + '\n' + importmap + '<script type="module">\n' + js + '\n</script>\n'
out = os.path.join(root, 'site', 'index.html')
open(out, 'w', encoding='utf-8').write(html)
print(out, len(html) // 1000, 'KB')
