import os
import shutil
import subprocess


def backup_original_files(files):
    for f in files:
        shutil.copy(f, f + ".bck")


def replace_in_file(file_name, regexes, suffix=''):
    for regex in regexes:
        subprocess.check_call(['sed', '-i' + suffix, regex, file_name])


def restore_file(file_name, suffix='.bck'):
    backup = file_name + '.bck'
    if (os.path.exists(backup)):
        shutil.move(backup, file_name)


def check_env():
    if not 'CUBROOT' in os.environ:
        raise EnvironmentError(
              "CUB (https://github.com/NVlabs/cub) should be downloaded and "
              "the environment variable 'CUBROOT' set to where it lies. "
              "No installation required: it's header only.")
    try:
        subprocess.check_output(["nvcc", "--version"])
    except:
       raise EnvironmentError("nvcc should be available on the system PATH")


def setup_files(makefile, config_h):
    replace_in_file(makefile,
                    [r's/^\s*KNC_BUILD\s*/# &/',
                     r's/^\s*ICC.*/ICC := icpc/',
                     r's/^#\(USE_CUDA := yes\)/\1/'])
            
    replace_in_file(config_h,
                    [r's/\(\s*constexpr\s*int\s*nEtaPart\s*= \).*\(\s*;.*\)/\1 1\2/',
                     r's/\(\s*constexpr\s*int\s*maxCandsPerSeed\s*= \).*\(\s*;\s*\/\/.*\)/\1 8\2/'])


def cleanup(makefile, config_h):
    restore_file(makefile)
    restore_file(config_h)


def make_distclean(rootdir):
    saved_path = os.path.abspath(os.path.curdir)
    os.chdir(rootdir)
    subprocess.check_call(["make", "-f", "Makefile", "distclean"])
    os.chdir(saved_path)


def make_all(rootdir):
    saved_path = os.path.abspath(os.path.curdir)
    os.chdir(rootdir)
    subprocess.check_call(["make", "-j", "24", "-f", "Makefile", "all"])
    os.chdir(saved_path)


def get_file_path_safe(rootdir, file_name):
    path = os.path.join(rootdir, file_name)
    if not os.path.isfile(path):
        raise RuntimeError('%s should exist, check the path' % path)
    return path


def get_file_paths(rootdir):
    makefile = get_file_path_safe(rootdir, 'Makefile.config')
    config_h = get_file_path_safe(rootdir, 'Config.h')
    return (makefile, config_h)


