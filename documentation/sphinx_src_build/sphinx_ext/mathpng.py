#!/usr/bin/python
import os
try:
    from hashlib import md5
except ImportError:
    from md5 import md5
    
from docutils import nodes
from docutils.writers.html4css1 import HTMLTranslator
from sphinx.latexwriter import LaTeXTranslator

# Define LaTeX math node:
class latex_math(nodes.General, nodes.Element):
    pass

def math_role(role, rawtext, text, lineno, inliner,
              options={}, content=[]):
    i = rawtext.find('`')
    latex = rawtext[i+1:-1]
    node = latex_math(rawtext)
    node['latex'] = latex
    return [node], []


try:
    from docutils.parsers.rst import Directive
except ImportError:
    # Register directive the old way:
    from docutils.parsers.rst.directives import _directives
    def math_directive(name, arguments, options, content, lineno,
                       content_offset, block_text, state, state_machine):
        latex = ''.join(content)
        node = latex_math(block_text)
        node['latex'] = latex
        return [node]
    math_directive.arguments = None
    math_directive.options = {}
    math_directive.content = 1
    _directives['math'] = math_directive
else:
    class math_directive(Directive):
        has_content = True
        def run(self): 
            latex = ' '.join(self.content)
            node = latex_math(self.block_text)
            node['latex'] = latex
            return [node]
    from docutils.parsers.rst import directives
    directives.register_directive('math', math_directive)

def setup(app):
    app.add_node(latex_math)
    app.add_role('math', math_role)

    # Add visit/depart methods to HTML-Translator:
    def visit_latex_math_html(self, node):
        source = self.document.attributes['source']
        self.body.append(latex2html(node, source))
    def depart_latex_math_html(self, node):
            pass
    HTMLTranslator.visit_latex_math = visit_latex_math_html
    HTMLTranslator.depart_latex_math = depart_latex_math_html

    # Add visit/depart methods to LaTeX-Translator:
    def visit_latex_math_latex(self, node):
        inline = isinstance(node.parent, nodes.TextElement)
        if inline:
            self.body.append('$%s$' % node['latex'])
        else:
            self.body.extend(['\\begin{align}\\begin{split}',
                              node['latex'],
                              '\\end{split}\\end{align}'])
    def depart_latex_math_latex(self, node):
            pass
    LaTeXTranslator.visit_latex_math = visit_latex_math_latex
    LaTeXTranslator.depart_latex_math = depart_latex_math_latex

from os.path import isfile,isdir
# LaTeX to HTML translation stuff:
def latex2html(node, source):
    inline = isinstance(node.parent, nodes.TextElement)
    latex = node['latex']
    print latex
    name = 'math-' + md5(latex).hexdigest()[-10:]
    # Always called from the document root directory.
    # Place the math images in the mathimg subdirectory 
    # of the document root directory.
    print os.getcwd
    path=os.path.join(os.getcwd(),'build/html/','mathimg')
    print 'path:',path
    if not isdir(path):
        os.mkdir(path)
    if not isfile('%s/%s.png' % (path,name)):
        f = open('math.tex', 'w')
        f.write(r"""\documentclass[12pt]{article}
                    \pagestyle{empty}""")
        g=open('source/math-defs.tex','r')
        f.write(g.read())
        g.close()
        f.write(r"\begin{document}")
        if inline:
            f.write('$%s$' % latex)
        else:
            f.write(r'\begin{align*} %s \end{align*}' % latex)
        f.write('\end{document}')
        f.close()
        os.system('latex --interaction=nonstopmode math.tex > /dev/null')
        os.system('dvipng -bgTransparent -Ttight --noghostscript -l10 ' +
                  '-o %s/%s.png math.dvi > /dev/null' % (path,name))
    if inline and '_' in latex:
        align = 'align="absmiddle" '
    else:
        align = ''
    if inline:
        cls = ''
    else:
        cls = 'class="center" '
    print 'source',source
    path=source.split('/source/')[-1].count('/')*'../'+'mathimg'
    return '<img src="%s/%s.png" %s%s/>' % (path, name, align, cls)
