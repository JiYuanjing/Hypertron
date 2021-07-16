 #!/bin/bash
 Address=$PWD
 echo $Address

filelist="mixeventTree.list"
mkdir log err csh report out_KF out_ME out_SE out_RT list
star-submit-template -template minitree_mixevent.xml -entities listOfFiles=${filelist},dir=$Address
filelist="rotateTree.list"
star-submit-template -template minitree_rotate.xml -entities listOfFiles=${filelist},dir=$Address
