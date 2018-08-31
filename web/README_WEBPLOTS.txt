##########################
## IMPORTANT DISCLAIMER ##
##########################

Please read the README.md before launching any validation scripts! It has useful info on scripts within this directory.

===================================

################################################
## Setting up a website with an EOS userspace ##
################################################

## Intro ##

First, ensure you have an LXPLUS account. If so, upon a standard login into LXPLUS will take you to your /afs home directory: /afs/cern.ch/user/${FIRST_LETTER}/${USERNAME}, where the first letter of your username is the variable ${FIRST_LETTER}, e.g. Kevin's ${USERNAME} is "kmcdermo", so ${FIRST_LETTER} is "k". Your corresponding /eos user space is: /eos/user/${FIRST_LETTER}/${USERNAME}, which is allotted 1TB of storage. If you are a member of CMS, you may have additional storage on /eos at LXPLUS or the FNAL LPC. However, for the purposes of building a website, use the /eos/user space for storing files for the web. It should be noted that your /eos/user space is also linked to CERNBox, and files can be shared quickly through this tool. Additional info on EOS and how to work with it can be found here: https://cern.service-now.com/service-portal/article.do?n=KB0001998

Now, you have to determine if you are either A) looking to migrate your personal website at CERN from AFS to EOS, or B) create a new personal website on EOS. In order to do option B), you cannot already have a personal website with CERN. So by default, if you want to setup an website with an EOS space at CERN, and you already have an AFS website, you will need to choose option A).

## Follow the procedure for setting up a website ##

1) Make sure to sign into your CERN account and go to this help page: https://cern.service-now.com/service-portal/article.do?n=KB0004096
2) The link listed as "Please click on the link" does not actually work, but the images do. Download and follow the instructions listed in eosuser-web-[1-4].png.
3) At this point, you will need to request a website from CERN. Follow the branches A or B below

   A) Migrating to EOS from AFS
      a) First, have a look at this link and watch the video: https://cds.cern.ch/record/2286177?ln=en , or have a look at this document: https://indico.cern.ch/event/661564/attachments/1512296/2358778/Migrate_Website_from_AFS_to_EOS.DOCX
      b) Follow the instructions in the video, ensuring to read the text once you click on the button: "Migrate to EOS". You will have to copy your old files over to /eos/user if you want the transition to be seamless.
   B) Brand new personal site at CERN
      a) Follow the instructions listed in eosuser-web-10.png.

4) While waiting for the request, you will need to setup directory browsing and persimissions: https://espace.cern.ch/webservices-help/websitemanagement/ConfiguringAFSSites/Pages/default.aspx
   A) If you do not care who views your website, at a mimumum, you need to follow the instructions on "Directory Browsing for AFS / EOS sites"
      a) cd /eos/user/${FIRST_LETTER}/${USERNAME}/${dir}, where ${dir} == "www" (but does not have to be, just whatever you named the directory that you requested for a website)
      b) touch .htaccess
      c) [open .htaccess and paste the following] "Options +Indexes" [without quotes]
   B) If you want to set restrictions on who can view your website, you can modify your .htaccess file by following along in "Access control for AFS / EOS web sites"

5) Move/copy new plots to /eos/user/${FIRST_LETTER}/${USERNAME}/${dir}, and use ./copyphp.sh ${subdir} to copy index.php into every directory below ${subdir}.

===================================

#################################
## WARNING: DEPRECATED BY CERN ##
#################################

A) Request CERN website from websites.cern.ch
   --> set website to point to AFS directory
   --> make website match username
   --> set address to /afs/cern.ch/user/${FIRST_LETTER}/${USERNAME}/${dir}
   --> make sure ${dir} exists!

B) While waiting for request, do the follow commands in one directory above ${dir}

  1) fs setacl ${dir} webserver:afs read
  2) afind ${dir} -t d -e "fs setacl -dir {} -acl webserver:afs read"
  3) cd ${dir}
  4) touch .htaccess
  5) [open .htaccess and paste the following] "Options +Indexes" [without quotes]

C) Then copy in really the very useful index.php into ${dir}

D) Once set up and website is live, copy plots and directories into {dir}
E) cd ${dir}
F) ./makereadable.sh ${subdir}, for every subdir. If launched from the top-level directory ${subdir}, it will handle the searching of subdirs.

--> As an aside, there are two other directories on LXPLUS every user has access to: /afs/cern.ch/ubackup/${FIRST_LETTER}/${USERNAME} and /afs/cern.ch/work/${FIRST_LETTER}/${USERNAME}. "ubackup" is a backup of 24h snapshots of "user", while "work" is not backed up but users can request up to 100 GB of space. The max for "user" directories is 10 GB upon request.
