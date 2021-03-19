

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('pyrad', parent_package, top_path)
    config.add_subpackage('util')
    config.add_subpackage('io')
    config.add_subpackage('proc')
    config.add_subpackage('prod')
    config.add_subpackage('graph')
    config.add_subpackage('flow')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
