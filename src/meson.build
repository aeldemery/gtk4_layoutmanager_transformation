gnome = import('gnome')
math_dep = math_project.get_variable ('math_dep')

gtk4_layoutmanager_transformation_sources = [
  'Main.vala',
  'MainWindow.vala',
  'GlobeWidget.vala',
  'GlobeTransformationLayout.vala',
]

gtk4_layoutmanager_transformation_deps = [
  dependency('gio-2.0', version: '>= 2.50'),
  dependency('gtk4', version: '>= 3.99.1'),
  dependency('gsl'),
  meson.get_compiler ('c').find_library ('m', required: false),
  meson.get_compiler ('vala').find_library ('posix'),
  meson.get_compiler ('vala').find_library ('MathHelper', dirs: math_vapi),
  math_dep,
]

gtk4_layoutmanager_transformation_sources += gnome.compile_resources('gtk4_layoutmanager_transformation-resources',
  'gtk4_layoutmanager_transformation.gresource.xml',
  c_name: 'gtk4_layoutmanager_transformation'
)

executable(
    'gtk4_layoutmanager_transformation',
    gtk4_layoutmanager_transformation_sources,
    vala_args: ['--vapidir', 'subprojects' / 'MathHelper'],
    dependencies: gtk4_layoutmanager_transformation_deps,
    install: true,
)
