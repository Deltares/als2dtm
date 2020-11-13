using Pkg
@info "Downloading and installing all required packages..."
Pkg.instantiate()
@info "Precompiling the als2dtm package..."
using als2dtm
