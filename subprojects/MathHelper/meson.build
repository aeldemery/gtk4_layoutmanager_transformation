project ('MathHelper', 'c', version: '0.2')

lib_math_sources = files (
    'four_point_transform.c',
    'singular_value_decomposition.c',
)

lib_math_include = include_directories ('.')

lib_math_dependencies = [
    dependency ('graphene-gobject-1.0'),
    dependency ('glib-2.0'),
    meson.get_compiler ('c').find_library ('m', required: false),
]

lib_math = static_library (
    'MathHelper', 
    sources: lib_math_sources,
    dependencies: lib_math_dependencies,
)

math_dep = declare_dependency (include_directories: lib_math_include, link_with: lib_math)
