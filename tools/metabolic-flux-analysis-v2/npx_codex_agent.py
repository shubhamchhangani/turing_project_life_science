"""Launch Harbor's isolated Codex evaluator through the requested npx CLI."""
from harbor.agents.installed.codex import Codex


class NpxCodex(Codex):
    async def exec_as_agent(self, environment, command, **kwargs):
        if "codex exec " in command:
            # Harbor installs this pinned package before the trial. npx selects
            # that same installed binary; no host task directory is mounted.
            command = command.replace(
                "codex exec ", "npx --yes @openai/codex@0.164.0 exec ", 1
            )
        return await super().exec_as_agent(environment, command, **kwargs)
