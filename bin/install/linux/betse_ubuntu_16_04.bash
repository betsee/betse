#!/usr/bin/env bash
# ====================[ betse_ubuntu_16_04.bash           ]====================
#
# --------------------( SYNOPSIS                           )--------------------
# Bash shell script automating the installation of BETSE and all mandatory and
# optional dependencies thereof (e.g., NumPy, SciPy) for all Ubuntu Linux
# releases newer than or equal to Ubuntu 16.04 (Xenial Xerus).
#
# --------------------( USAGE                              )--------------------
# ./betse_ubuntu_16_04.bash [INSTALL_DIR]
#
# Installs both BETSE and BETSE to the passed installation directory,
# defaulting to the common "${HOME}/Applications" directory if unpassed. This
# directory and all parent directories of this directory will be implicitly
# created as needed (i.e., if they do *NOT* currently exist).
#
# --------------------( DOWNLOAD                           )--------------------
# While this script may be run directly from a cloned Git repository, doing so
# would defeat the utility of this script -- which automates the cloning of this
# repository in a predictably structured manner, among other useful installation
# tasks. Instead, users are recommended to directly download this script from
# BETSE's GitLab-hosted project to the local filesystem and run that.
#
# Since all Ubuntu releases come preinstalled with the "wget" command, the
# following "one-liner" suffices to install BETSE from any open terminal:
#
#     wget https://gitlab.com/betse/betse/raw/master/bin/install/linux/betse_ubuntu_16_04.bash && bash betse_ubuntu_16_04.bash

#FIXME: Adequately test this installation script and, if seemingly working,
#revise our "README.rst" instructions accordingly.

#FIXME: *UGH.* The scipy.misc.imread() has been officially deprecated and will
#be removed as of SciPy 1.2.0. SciPy now recommends the "imageio" package as a
#replacement dependency providing an equivalent function. We therefore need to
#perform the following:
#
#* Define a new BETSE-specific "bin/install/conda.bash" shell script
#  performing an Anaconda-based installation -- complete with Git repository
#  cloning stripped straight out of our existing Ubuntu installer. This
#  script should probably *NOT* attempt to install Anaconda but simply fail
#  with a fatal error in the absence of "conda" in the current ${PATH}.
#  Note that this script should obviate the need for users to manually
#  download and install Anaconda now via automated commands resembling:
#
#  # This isn't the entire story, of course. This shell script also needs to
#  # re-source "~/.bashrc" due to Anaconda changes. See this (for details):
#  #     http://web.workstate.com/blog/workstate-codes-openai-gym-on-windows-ubuntu-bash-shell-yes-it-works
#  $ wget https://repo.continuum.io/archive/Anaconda3-4.2.0-Linux-x86_64.sh
#  $ bash Anaconda3-4.2.0-Linux-x86_64.sh
#* Document both our Ubuntu >= 16.04 and Anaconda installers (in that order,
#  as Ubuntu users would certainly prefer a process integrating with their
#  existing package manager) at the head of "README.md".
#* Document more explicitly that Windows users should *ABSOLUTELY* install
#  Bash-on-Ubuntu-on-Windows. This should literally be the first item in
#  this Markdown list of installation instructions, as it then permits
#  Windows users to install BETSE via either:
#  * The Ubuntu shell script (strongly recommended).
#  * The Anaconda shell script (feasible but less recommended).
#  Actually, we need to create a Windows 10-specific installation script.
#  Attempting to leveraging an existing non-Windows installation script to
#  install under Windows is a recipe destined for obvious failure. Still,
#  this script should strongly resemble our Anaconda installation script,
#  but will at least need to be a superset of this script (e.g., due to the
#  need to install Xming or a similar display server). Ideally, our
#  Windows 10-specific installation script would implicitly download and
#  source our Anaconda installation script -- but perhaps there will be too
#  many differences between the two to support sourcing? Let's find out,
#  shall we. Further details on Anaconda under Ubuntu under Windows here:
#      http://web.workstate.com/blog/workstate-codes-openai-gym-on-windows-ubuntu-bash-shell-yes-it-works
#  Note lastly that the prior link leverages Xming. That doesn't work, in
#  our case, as Xming is proprietary non-free software. The sane alternative
#  is VcXsrv. Ideally, our Windows installation script will implicitly
#  download and install VcXsrv as well. Is that even feasible? *sigh*
#  VcXsrv lives at:
#      https://sourceforge.net/projects/vcxsrv
#  Ah-ha! Found it. The following Bash shell script automates installation
#  of VcXsrv in Ubuntu on Windows. Pretty slick, honestly:
#      https://github.com/patrick330602/wslu/blob/master/wslpkg#L83
#  Also, the following installation snippet will probably prove useful:
#      # If "~/.bashrc" does *NOT* already contain a line exporting an X.org
#      # display socket, append a line exporting the default socket:
#      grep '^export DISPLAY' ~/.bashrc ||
#          echo 'export DISPLAY=localhost:0.0' >> ~/.bashrc
#  Note lastly that our Windows installation script will probably need to
#  explicitly install an X.org package via "apt" -- say, "x11-dev", maybe?

# ....................{ BASH                               }....................
# Enable strict mode, terminating the script with non-zero exit status on the
# first command or pipeline failing with non-zero exit status.
set -o errexit

# ....................{ FUNCTIONS ~ print                  }....................
# note(info_message: str) -> None
#
# Print the passed informational message to standard output.
note() {
    echo "${BASH_SOURCE[1]}:" "${@}"
}


# info(info_message: str) -> None
#
# Print the passed informational message to standard output prefixed by a
# newline.
info() {
    echo
    note "${@}"
}


# die(error_message: str) -> None
#
# Print the passed error message to standard error prefixed by a newline and
# terminate the current shell script with the usual failure exit status.
die() {
    echo "${BASH_SOURCE[1]}:" "${@}" >&2
    exit 1
}

# ....................{ FUNCTIONS ~ install                }....................
# install_apt_package(*package_names: str) -> None
#
# Install all Ubuntu packages with the passed names via the "apt" command.
#
# For convenience, this function additionally upgrades these packages and all
# transitive dependencies of these packages if already installed.
install_apt_package() {
    (( ${#} >= 1 )) || die 'One or more arguments expected.'

    note "Installing Ubuntu package(s) \"${*}\"..."
    sudo apt install --yes "${@}"
}


# install_pip_package(*package_names: str) -> None
#
# Install all Python packages with the passed names via the "pip3" command from
# the remote PyPI-hosted distributions of the same names.
#
# For convenience, this function additionally upgrades these packages and all
# transitive dependencies of these packages if already installed.
install_pip_package() {
    (( ${#} >= 1 )) || die 'One or more arguments expected.'

    note "Installing Python package(s) \"${*}\"..."
    sudo --set-home pip3 --disable-pip-version-check install --upgrade "${@}"
}


# install_pip_git_repo(repo_url: str, repo_dirname: str) -> None
#
# Install the Python project hosted at the remote Git repository with the passed
# URL to the local directory with the passed absolute or relative path in an
# editable manner, preventing repository synchronization issues.
install_pip_git_repo() {
    (( ${#} == 2 )) || die 'Two arguments expected.'

    local repo_url="${1}"
    local repo_dirname="${2}"

    # If this directory already exists...
    if [[ -d "${repo_dirname}" ]]; then
        # If this directory is *NOT* a Git repository, fail.
        [[ -d "${repo_dirname}/.git" ]] ||
            die "Directory \"${repo_dirname}\" not a Git repository."

        # Synchronize this local repository against remote changes.
        note 'Updating Git repository...'
        pushd "${repo_dirname}"
        git pull
        popd
    # Else, clone this remote repository to this local directory.
    else
        note 'Cloning Git repository...'
        git clone "${repo_url}" "${repo_dirname}"
    fi

    # If this directory contains no "setup.py" script, fail.
    [[ -f "${repo_dirname}/setup.py" ]] ||
        die "\"${repo_dirname}\" not a Python project."

    # Install this Python project editably.
    note 'Installing Python project...'
    install_pip_package --editable "${repo_dirname}"
    # pushd "${repo_dirname}"
    # sudo python3 setup.py develop
    # popd
}

# ....................{ FUNCTIONS ~ os                     }....................
# is_os_linux_ubuntu() -> int
#
# Success only if the current platform is Ubuntu Linux.
is_os_linux_ubuntu() {
    python - << '__HEREDOC'
import platform, sys

# If the Linux-specific platform.linux_distribution() function is available
# *AND* the first item of the 3-tuple "(distname, version, id)" returned by
# this function is "Ubuntu", the current platform is Ubuntu Linux. In this
# case, the active Python process is terminated with the success exit code.
# Since the platform.system() function returns the low-level kernel name
# "Linux", this function is *NOT* called here.
if (hasattr(platform, 'linux_distribution') and
    platform.linux_distribution()[0] == 'Ubuntu'):
    sys.exit(0)
# Else, the current platform is *NOT* Ubuntu Linux. In this case, the active
# Python process is terminated with a failure exit code.
else:
    sys.exit(1)
__HEREDOC
}

# ....................{ PREAMBLE                           }....................
# Print a quasi-informative preamble.
note 'Welcome to the BETSE installer for Ubuntu Linux 16.04 and newer.'
echo

# If sudo privelages have already expired for the current user, (re)cache these
# privelages in a human-readable manner. The default prompt printed by the
# "sudo" command (e.g., "[sudo] password for ${USERNAME}:") is arguably
# non-human-readable for those unfamiliar with POSIX shell environments.
if ! sudo -S true </dev/null 2> /dev/null; then
    note 'Please enter your user password.'
    sudo true
fi

# ....................{ CHECKS                             }....................
# If the current platform is *NOT* Ubuntu Linux, fail.
note 'Detecting current platform...'
is_os_linux_ubuntu || die 'Current platform not Ubuntu Linux.'

# ....................{ ARGS                               }....................
note 'Parsing script arguments...'

# Absolute or relative path of the parent directory containing the directories
# to which BETSE and BETSE will be installed, defaulting to a general-purpose
# directory if unpassed.
install_dirname="${1-${HOME}/Applications}"

# Absolute or relative path of the directory to which BETSE will be installed.
betse_dirname="${install_dirname}/betse"

# Absolute or relative path of the directory to which BETSE will be installed.
betse_dirname="${install_dirname}/betse"

# Create this directory and all parent directories of this directory as needed.
note "Installing BETSE to: ${install_dirname}"
mkdir --parents "${install_dirname}"

# ....................{ DEPENDENCIES                       }....................
# Install subsequently required core commands (e.g., "git", "pip3").
note 'Installing project and package managers...'
install_apt_package git python3-pip

# ....................{ DEPENDENCIES ~ betse               }....................
# Note that the entirety of this section is derived from the Debian-specific
# installation instructions for BETSE at:
#     https://gitlab.com/betse/betse/blob/master/doc/md/INSTALL.md#debian

# Install all mandatory BETSE-specific dependencies.
info 'Installing mandatory BETSE dependencies...'
install_apt_package \
    python3-dev python3-dill python3-matplotlib \
    python3-numpy python3-pil python3-pip python3-scipy python3-setuptools \
    python3-six tcl tk

# The living "ruamel.yaml" is strongly preferred to the dead PyYAML.
install_pip_package ruamel.yaml

# Install an OpenBLAS-optimized scientific stack. While the SourceForge-hosted
# ATLAS project would also apply, the GitHub-hosted OpenBLAS project is
# (unsurprisingly) significantly better maintained.
info 'Installing OpenBLAS-optimized scientific stack...'
install_apt_package build-essential libopenblas-dev
sudo update-alternatives --set libblas.so.3 /usr/lib/openblas-base/libblas.so.3
sudo update-alternatives --set liblapack.so.3 /usr/lib/lapack/liblapack.so.3

#FIXME: Add installation instructions to "INSTALL.md" documenting installation
#of FFMpeg and LibAV on all supported platforms.

# Install all optional BETSE-specific dependencies.
#
# Note that:
#
# * The "libav-tools" package provides the "avconv" command required for
#   exporting compressed videos of simulation phase runs.
# * The "python3-networkx" package provides NetworkX 1.11, which breaks
#   backward compatibility with respect to PyDot support required by BETSE.
#   Hence, the next-most-recent NetworkX release is installed manually.
info 'Installing optional BETSE dependencies...'
install_apt_package graphviz libav-tools python3-pydot
install_pip_package 'networkx==1.10'

# ....................{ BETSE                              }....................
# Install the live version of BETSE in an editable manner. As of this writing,
# most BETSE users are currently also developers and hence prefer the
# developer-oriented installation performed here.
info 'Installing BETSE...'
install_pip_git_repo 'https://gitlab.com/betse/betse.git' "${betse_dirname}"

# ....................{ INSTRUCTIONS                       }....................
# For usability, provide rudimentary usage instructions.
info "Congratulations! BETSE has been successfully installed.'

BETSE is now runnable at the command line with the following command:
    betse
"
