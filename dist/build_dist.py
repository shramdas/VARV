#!/usr/bin/env python
import os
import os.path as path
import sys

# Move up to the runepacts root.
this_dir = path.dirname(path.realpath(sys.argv[0]))
root_dir = path.join(this_dir,"..")

print "Changing into: %s" % root_dir
os.chdir(root_dir)
sys.path.append(root_dir)

from runepacts.main import PROG_VERSION

tarfile = "runepacts_%s.tgz" % PROG_VERSION;

if os.path.exists(tarfile):
  print "Removing existing tar.."
  os.remove(tarfile);

command = "tar zcfh %s" % tarfile;
paths = """
bin/runepacts
bin/setup.py
bin/virtualenv-runepacts.py
bin/virtualenv_support
runepacts
README.md
""".split();

args = [];

print "Options enabled: "

if 1:
  print ".. adding runepacts/ to beginning of each file path";
  args.append("--xform 's|^|runepacts/|'");

if 1:
  print ".. removing git";
  args.append("--exclude *git*");

if 1:
  print ".. removing pycharm files";
  args.append("--exclude *.idea*");

if 1:
  print ".. removing .pyo and .pyc files";
  args.append("--exclude *.pyo*");
  args.append("--exclude *.pyc*");

if 1:
  print ".. removing GoT2D files";
  args.append("--exclude *GoT2D*");

if 1:
  print ".. removing FUSION GWAS catalog files";
  args.append("--exclude *gwascat_fusion*")

final = "%s %s" % (command," ".join(paths + args));

print "Packing.."
os.system(final);

