# How do I choose the right tutorial for my platform?

There is a [platform tutorial](../vignettes/intro.md#platform-tutorials) for each of these platforms: 10x Xenium, Vizgen MERSCOPE, CosMx SMI, 10x Visium HD, Stereo-seq, Seq-Scope, and Pixel-seq, plus a [Generic tutorial](../vignettes/subregion_tutorials/generic.md) for any other platform's transcript table. Each starts from a small dataset, shows the single `run_together` command, and then walks through the individual commands.

For your own data, see the platform's page under [Platforms](../reference/supported_platforms.md): it lists the exact files `run_together` looks for and example commands. Built-in profiles exist for 10x Xenium, 10x Visium HD, CosMx SMI, MERSCOPE (MERFISH), Illumina, Stereo-seq, and Seq-Scope, plus a `generic` profile for anything else.

If your platform is not listed above, see:

- [My platform is not in the supported list. Can I still use `CartLoader`?](./unsupported_platform.md)
