# V2 frontier trial status — DO NOT SUBMIT

The fresh `npx @openai/codex@0.164.0` trial using the existing ChatGPT login
completed successfully on GPT-6-sol, max effort. The agent saw only the task
instructions and agent-visible input files in an isolated Docker container.
It independently produced both artifacts and received reward **1**. Pytest
reported **16 passed**. The agent's full-data weighted SSR was 1961.0084,
below the allowed 2110. Agent execution took 11 minutes 18 seconds.

This candidate therefore fails the user's difficulty requirement. No new ZIP
was created after this result. The older candidate ZIP is not recommended for
submission. Original logs and artifacts remain under
`jobs/mfa-v2-frontier/gpt6-sol-max-v2-npx-chatgpt`.

The static preflight passed 37 checks. The offline artifact regression suite
accepted four valid cases and rejected all 25 invalid cases, including the
previously reported PGI exchange loophole, decimal-field whitespace, Unicode
digits, malformed headers, and FIFO files. Those checks establish package and
verifier behavior; they do not establish frontier-model difficulty.

## Earlier API attempt

The authorized GPT-6-sol max-effort run ended with `ApiUsageLimitError`:
“Quota exceeded. Check your plan and billing details.”

The agent log contains a turn start followed by the quota error, with no
substantive solution attempt. Harbor recorded reward 0 because no valid
artifacts were produced. This is an infrastructure failure and must not be
counted as evidence that the task defeated the model. Token usage and cost were
not reported; no claim about the actual charge is made.

The v2 task's offline reference solution passed all 16 verifier tests. At that
point its frontier difficulty was unmeasured; the later npx trial above resolved
that uncertainty for GPT-6-sol on the local Codex harness.

The API-key run was not retried. The subsequent user-requested npx run used
the existing ChatGPT login instead; it did not reuse the API-key route.

See `jobs/mfa-v2-frontier/gpt6-sol-max-v2/evaluation-status.json` for the
classification and pointer to the unmodified raw result.
