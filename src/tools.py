import os
import shutil


def tool(name):
    """Command used to call an external program.

    A copy sitting in the working directory wins, which is how the standalone
    and Docker distributions ship their binaries, so they behave as before.
    Otherwise the program is looked up on PATH (conda, a module system, a
    system install), which lets MeTarEnz run from any directory.
    """
    if os.path.exists(name):
        return './' + name
    return shutil.which(name) or name
