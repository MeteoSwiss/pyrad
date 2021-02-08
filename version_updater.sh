
#!/bin/bash
echo Do you want to release a new pyART version y/n?
read updatepyart
if [[ $updatepyart == "y" ]]; then
echo "Please provide version number, for example 0(major).5(minor).0(micro)"
echo "What is the major version?"
read major
echo "What is the minor version?"
read minor
echo "What is the micro version?"
read micro

echo "Version number is $major.$minor.$micro"
echo "Updating setup.py"
sed -i "/MAJOR =/c\MAJOR = $major" ./src/pyart/setup.py
sed -i "/MINOR =/c\MINOR = $minor" ./src/pyart/setup.py
sed -i "/MICRO =/c\MICRO = $micro" ./src/pyart/setup.py

echo "Updating /doc/source/conf.py"
sed -i "/version = u/c\version = u'$major.$minor'" ./src/pyart/doc/source/conf.py
sed -i "/release = u/c\release = u'$major.$minor.$micro'" ./src/pyart/doc/source/conf.py

echo "Done! You can now create and push the new release with:"
echo "cd ./src/pyart/" 
echo "git tag v$major.$minor"
echo "git push origin v$major.$minor."
fi

echo Do you want to release a new pyRAD version y/n?
read updatepyrad
if [[  "$updatepyrad" == "y" ]]; then
echo "Please provide version number, for example 0(major).5(minor).0(micro)"
echo "What is the major version?"
read major
echo "What is the minor version?"
read minor
echo "What is the micro version?"
read micro

echo "Version number is $major.$minor.$micro"
echo "Updating setup.py"
sed -i "/MAJOR =/c\MAJOR = $major" ./src/pyrad_proc/setup.py
sed -i "/MINOR =/c\MINOR = $minor" ./src/pyrad_proc/setup.py
sed -i "/MICRO =/c\MICRO = $micro" ./src/pyrad_proc/setup.py

echo "Updating /doc/source/conf.py"
sed -i "/version = u/c\version = u'$major.$minor'" ./doc/source/conf.py
sed -i "/release = u/c\release = u'$major.$minor.$micro'" ./doc/source/conf.py

echo "Done! You can now create and push the new release with:"
echo "git tag v$major.$minor"
echo "git push origin v$major.$minor."
fi
