# -*- coding: utf-8 -*-
"""
    sphinxcontrib.dir
    ~~~~~~~~~~~~~~~~~~
    :copyright: Copyright 2012 by Takeshi KOMIYA
    :license: BSD, see LICENSE for details.
"""

import os

from docutils import nodes, utils
from docutils.parsers.rst import directives

try:
    import sphinx.util.compat
except ImportError:
    import sys
    import types
    import sphinx.util
    import docutils.parsers.rst
    class compat(types.ModuleType):
        Directive = docutils.parsers.rst.Directive
    sphinx.util.compat = compat('sphinx.util.compat')
    sys.modules['sphinx.util.compat'] = sphinx.util.compat

Directive = sphinx.util.compat.Directive


class directory(nodes.General, nodes.Element):
    pass


class DirDirective(Directive):
    """Directive for embedding google-maps"""

    has_content = False
    required_arguments = 0
    optional_arguments = 1
    final_argument_whitespace = True
    option_spec = {
    }

    def run(self):
        node = directory()
        prefix = os.path.abspath(os.path.dirname(__file__)+'/../../../nullspace_optimizer/')
        node['files'] = [self.arguments[0]+'/'+x for x in os.listdir(prefix+'/'+self.arguments[0]) if x.endswith('.py')]
        node['files'].sort()

        document = self.state.document
        if 'files' not in node:
            msg = 'dir directive needs argument for searching files'
            return [document.reporter.warning(msg, line=self.lineno)]

        return [node]


def visit_directory_node(self, node):
    self.body.append('<ul>')

    prefix_url = "https://gitlab.com/florian.feppon/null-space-optimizer/-/blob/public-master/nullspace_optimizer/"  
    for file in node['files']:
        icon = r"""<svg aria-hidden="true" class="external-link-icon" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg"><path d="M19 19H5V5h7V3H5a2 2 0 00-2 2v14a2 2 0 002 2h14c1.1 0 2-.9 2-2v-7h-2v7zM14 3v2h3.59l-9.83 9.83 1.41 1.41L19 6.41V10h2V3h-7z"></path></svg>"""
        self.body.append('<li><p><a class="reference external" href="'+prefix_url+'/'+file+'" rel="nofollow noopener">' 
                          + file +  icon + '</a></p></li>')

    self.body.append('</ul>')


def depart_directory_node(self, node):
    pass


def setup(app):
    app.add_node(directory,
                 html=(visit_directory_node, depart_directory_node))
    app.add_directive('dir', DirDirective)
