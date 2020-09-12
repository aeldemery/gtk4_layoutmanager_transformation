project('gtk4_layoutmanager_transformation', ['c', 'vala'], 
    version: '0.1.0',
    meson_version: '>= 0.53.0',
    default_options: [ 'warning_level=2',],
)

math_vapi = meson.current_source_dir() / 'subprojects' / 'MathHelper'

math_project = subproject ('MathHelper')

#add_project_arguments(['--vapidir', math_vapi], language: 'vala')

i18n = import('i18n')

subdir('data')
subdir('src')
subdir('po')

meson.add_install_script('build-aux/meson/postinstall.py')