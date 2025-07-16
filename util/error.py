import sys
import traceback
import pdb

def react(f):
    try:
        f()
    except:
        traceback.print_exc()
        pdb.post_mortem(sys.exc_info()[2])
