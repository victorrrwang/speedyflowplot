# Building speedyflowplot Package

This guide explains how to build the speedyflowplot package from the root folder.

## Building from Root Folder

### Option 1: Using R CMD build (command line)

From the root folder, run:

```bash
R CMD build pkg
```

This will create a `.tar.gz` file in the root folder. Then install it with:

```bash
R CMD INSTALL speedyflowplot_*.tar.gz
```

Or to build directly to the releases folder:

```bash
R CMD build pkg && mv speedyflowplot_*.tar.gz releases/
```

### Option 2: Using the Makefile

From the root folder, you can run:

```bash
cd pkg && make build
```

Or in one command:

```bash
(cd pkg && make build)
```

This will build the package and automatically move the `.tar.gz` file to the `releases/` directory.

### Option 3: Using R devtools (recommended for development)

From R, run:

```r
# From the root folder, you can use:
devtools::install("pkg")

# Or if you want to build without installing:
devtools::build("pkg")

# Build to a specific location
devtools::build("pkg", path = "releases")
```

### Option 4: Direct R CMD build with full path

```bash
R CMD build /path/to/speedyflowplot/pkg
```

## Quick Reference

- **Quick build**: `R CMD build pkg` (from root)
- **Development install**: `devtools::install("pkg")` (from R)
- **Using Makefile**: `cd pkg && make build` (automatically moves to releases/)

## Notes

The package includes C++ code that will be compiled during the build process. Make sure you have a C++ compiler configured for R (typically handled automatically on macOS with Xcode Command Line Tools).

## Additional Makefile Targets

If you're working in the `pkg` directory, the Makefile provides:

- `make build` - Build the package and move `.tar.gz` to `releases/`
- `make check` - Build and run R CMD check on the package in `releases/`
- `make install` - Build and install the package from `releases/`
- `make clean` - Remove build artifacts from `releases/` and compiled objects

