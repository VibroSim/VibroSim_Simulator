import subprocess
import os
import platform

def open_via_helper(filepath):
    """Open a file via the system viewer for the specified file. 
    Returns a handle to a cleanup function. You should call 
    handle() before exiting."""
    if hasattr(os,"startfile"):
        # On Windows the os module has a "startfile" function just for this
        os.startfile(filepath)
        return lambda: None  # startfile() fully detaches its subprocess
    elif platform.system() == "Darwin":
        # Mac OS X has a "open" command
        subproc = subprocess.Popen(["open",filepath])
        return subproc.communicate # caller should call this to wait for subproc to finish
    else:
        # Assume Linux/FreeDesktop
        subproc = subprocess.Popen(["xdg-open",filepath])
        return subproc.communicate # caller should call this to wait for subproc to finish
    pass
