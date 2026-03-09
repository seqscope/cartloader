# I am not sure whether to run with Docker or run locally. Which one should I choose?

Use this quick comparison to choose a run mode.

| Mode | Good | Bad | Best for |
|---|---|---|---|
| **Docker run** | Fastest setup; reproducible environment; fewer dependency/version issues; easier onboarding for new users | Docker install required; container startup overhead; file mounts/permissions can be confusing; less convenient for editing internals | First-time users, tutorials, shared team workflows, reproducible reruns |
| **Local run** | Direct access to local files/tools; easier debugging and profiling; no container boundary; easier to customize and extend | More setup effort; dependency conflicts are more likely; environments can drift across machines | Power users, developers, custom pipelines, tight integration with local tooling/HPC modules |

## Quick decision rule

- Choose **Docker** if your priority is getting started quickly and reproducibly.
- Choose **local** if your priority is deep customization, debugging, or development.

- Docker guide: [Run with Docker](../vignettes/quickstart/run_in_docker.md)
- Local guide: [Run Locally](../vignettes/quickstart/run_locally.md)
