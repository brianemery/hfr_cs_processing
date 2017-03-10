#hfr_cs_processing


''Main Git Commands''
{{{
# pull updates before making changes
git pull origin master

# pull updates from the origin's master repo, put your changes on top of master's
git pull --rebase origin master

# view the state of the repo
git status 

# stage changes (examples) (use this to group changes)
git add .
git add <directory>
git add <filename>

# make the commit
git commit -m "this is a message"

#Push commits to the central repository
git push origin master

# This will add deletes as well.
git add -u .


}}}

''Remove a directory and contents''
{{{
# Dont remove it from the command line, let git remove it
git rm -r the-directory
git commit -m "Remove duplicated directory"
}}}

''Set up a copy of the repo on a local machine''
{{{
git clone ssh://stokes/data/code_repo/m_files.git
}}}
