
#!/bin/bash
echo "This is a command line tool to help you release new pyart/pyrad version"
echo "Please confirm that you are on master branch and you have added and commit all modifications you want to include in the tag in the master branch!"
echo "Type enter to continue"
read continue

#####################################################
# PyART
#####################################################

dir="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
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

cd $dir

echo "Version number is $major.$minor.$micro"
echo "Updating setup.py"
sed -i "/MAJOR =/c\MAJOR = $major" $dir/src/pyart/setup.py
sed -i "/MINOR =/c\MINOR = $minor" $dir/src/pyart/setup.py
sed -i "/MICRO =/c\MICRO = $micro" $dir/src/pyart/setup.py

echo "Updating /doc/source/conf.py"
sed -i "/version = u/c\version = u'$major.$minor'" $dir/src/pyart/doc/source/conf.py
sed -i "/release = u/c\release = u'$major.$minor.$micro'" $dir/src/pyart/doc/source/conf.py

cd $dir/src/pyart/
echo "Done!"
echo "We will now push the version number updated on github"
git add --all
git commit -m "Changed version numbers"
git push origin master

echo "Do you want to push the new tag on github? y/n"
read pushtag
if [[ $pushtag == "y" ]]; then
git tag "v${major}.${minor}.${micro}"
git push origin "v${major}.${minor}.${micro}"
fi
echo "Done!"

echo "Do you want to create a new PyPI (pip) package? y/n"
read releasepypi
if [[ $releasepypi == "y" ]]; then
rm -r dist
python setup.py sdist
echo "Should it be a test package or a final (official) package? test/final"
read pkgtype
if [[ $pkgtype == "test" ]]; then
python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
echo "You can now install pkg with pip install --index-url https://test.pypi.org/simple/ --no-deps pyart-mch"
echo "Install a fresh virtual env to test it with"
echo "python –m venv pyart_testing"
echo "source pyart_testing/bin/activate"
fi
if [[ $pkgtype == "final" ]]; then
python -m twine upload dist/*
fi
fi
echo "Done!"
fi

#####################################################
# PyRAD
#####################################################

cd $dir
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

cd $dir

echo "Version number is $major.$minor.$micro"
echo "Updating setup.py"
sed -i "/MAJOR =/c\MAJOR = $major" ./src/pyrad_proc/setup.py
sed -i "/MINOR =/c\MINOR = $minor" ./src/pyrad_proc/setup.py
sed -i "/MICRO =/c\MICRO = $micro" ./src/pyrad_proc/setup.py

echo "Updating /doc/source/conf.py"
sed -i "/version = u/c\version = u'$major.$minor'" ./doc/source/conf.py
sed -i "/release = u/c\release = u'$major.$minor.$micro'" ./doc/source/conf.py

echo "We will now push the version number updated on github"
git add --all
git commit -m "Changed version numbers"
git push origin master

echo "Done!"
echo "Do you want to push the new tag on github? y/n"
read pushtag
if [[ $pushtag == "y" ]]; then
git tag "v${major}.${minor}.${micro}"
git push origin "v${major}.${minor}.${micro}"
fi

echo "Do you want to create a new PyPI (pip) package? y/n"
read releasepypi
if [[ $releasepypi == "y" ]]; then
cd $dir/src/pyrad_proc/
rm -r dist
python setup.py sdist
echo "Should it be a test package or a final (official) package? test/final"
read pkgtype
if [[ $pkgtype == "test" ]]; then
python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
echo "You can now install pkg with pip install --index-url https://test.pypi.org/simple/ --no-deps pyrad-mch"
echo "Install a fresh virtual env to test it with"
echo "python –m venv pyrad_testing"
echo "source pyrad_testing/bin/activate"
fi
if [[ $pkgtype == "final" ]]; then
python -m twine upload dist/*
fi
fi

fi

