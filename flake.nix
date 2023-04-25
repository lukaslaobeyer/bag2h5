{
  description = "bag2h5";

  outputs = { self, nixpkgs }: {
    devShells.x86_64-linux.default =
      let
        system = "x86_64-linux";
        pkgs = import nixpkgs { inherit system; };
      in
        with pkgs; mkShell {
          nativeBuildInputs = [ bazel_5 gdb ];
        };
  };

  nixConfig.bash-prompt = "\[nix-develop\]$ ";
}
