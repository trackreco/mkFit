To get CERN website up and running:

A) Request CERN website from websites.cern.ch
   --> set website to point to AFS directory
   --> make website match username
   --> set address to /afs/cern.ch/user/{u}/{username}/{dir}
   --> make sure {dir} exists!

B) While waiting for request, do the follow commands in one directory above {dir}

  1) fs setacl {dir} webserver:afs read
  2) afind {dir} -t d -e "fs setacl -dir {} -acl webserver:afs read"
  3) cd {dir}
  4) touch .htaccess
  5) [open .htaccess and paste the following] "Options +Indexes" [without quotes]

C) Then copy in really the very useful index.php into {dir}

D) Once set up and website is live, copy plots and directories into {dir}
E) cd {dir}
F) ./makereadable {subdir}, for every subdir
