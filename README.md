# SuperCompW1617
## Contents
  1.How to set up git on SuperMUC
  2.Bash file
## How to set up git on SuperMUC:

Only works for repositories of github.com, not for lrz gitlab repositories );

1. Log in to lxhalle like usual with:

  ```sh
  ssh <informatics-tum-login>@lxhalle.informatik.tu-muenchen.de
  ```

2.  Log in to SuperMUC with:
  ```sh
  ssh -R <port-number>:github.com:9418: <supermuc-login>@supermuc.lrz.de
  ```

  where port-number could be any number  between 10000 and 65535


3.  On supermuc:

  First load the git module with:
  ```sh
  module load git
  ```
  Then clone the git with
  ```sh
  git clone git://localhost:<port-number>/<git-repo>
  ```
  where port-number is the same as step 2 and git-repo is the name of the repo.
  For example if my repo in github has address "github.com/someone/myrepo.git" then git-repo would be someone/myrepo.git

  Just remember that in supermuc we can only pull from git but not push so you need to make the changes somewhere else.

## Bashfile

If you don't want to type the same commands every time you log in to SuperMUC then just create a file named **.bashrc** in the folder where you log in (you could do this with vim). There you should write the commands you want to be run every time you log in (one per line). This is done mainly for the **module load** or **module unload** commands for MPI and git.
