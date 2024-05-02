
#!/usr/bin/bash

if [ $1 == true ];
	then
	 rm -rf $2
	else
	echo "keeping work directory"
fi
