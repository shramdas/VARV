#!/usr/bin/env python

#===============================================================================
# Copyright (C) 2014 Shweta Ramdas, Ryan Welch, The University of Michigan
#
# runepacts is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# runepacts is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================

import virtualenv, textwrap

output = virtualenv.create_bootstrap_script(textwrap.dedent("""
import os, subprocess, shlex
def after_install(options, home_dir):
    etc = join(home_dir, 'etc')
    if not os.path.exists(etc):
        os.makedirs(etc)
    
    pip = join(home_dir,'bin','pip')
    
    cmd = "%s install -I pandas" % pip
    subprocess.call(shlex.split(cmd))
    
    cmd = "%s install -I numexpr" % pip
    subprocess.call(shlex.split(cmd))

"""))

f = open('virtualenv-varv.py', 'w').write(output)

