#controller = {'host':'128.32.142.141',  }

sshx = '~/.ssh/sshx'

#just autoconcat; nothing tuplish
ssh_options_string = ('-o UserKnownHostsFile=/dev/null '
                      '-o StrictHostKeyChecking=no '
                      '-i ~/.ssh/gsg-keypair '
                      '-l root')

push_kwargs = {'username':'root'}



