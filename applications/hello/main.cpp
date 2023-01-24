#include <iostream>
#include <forestclaw2d.h>
#include <fclaw2d_include_all.h>
#include <EllipticForestApp.hpp>

int main(int argc, char** argv) {
    std::cout << "Hello from fclaw-apps!" << std::endl;

    // Initialize ForestClaw
    fclaw_app_t* fclaw_app = fclaw_app_new(&argc, &argv, NULL);
    fclaw_global_essentialf("Hello from ForestClaw!");

    // Initialize EllipticForest
    EllipticForest::EllipticForestApp ef_app = EllipticForest::EllipticForestApp(&argc, &argv);
    ef_app.log("Hello from EllipticForest!");

    return EXIT_SUCCESS;
}