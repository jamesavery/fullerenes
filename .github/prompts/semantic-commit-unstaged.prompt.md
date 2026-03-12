---
name: Semantic Stage Unstaged Changes
description: "Inspect unstaged repository changes, group them into honest semantic commits, and stage the next coherent commit with a concise factual message proposal."
argument-hint: "Optional constraints for grouping or commit wording"
agent: "agent"
---
Review the repository's current unstaged changes and turn them into semantically related commit groups.

This repository does not allow the agent to create git commits directly. Instead, stage one coherent commit at a time and propose the commit message the user should use.

Follow this process:

1. Inspect the unstaged diff before changing the index.
2. Propose a commit plan that groups changed files and hunks by actual purpose, not by directory alone.
3. If the grouping is ambiguous, risky, or mixes unrelated work, stop and ask a short clarifying question before staging anything.
4. Otherwise, stage the next commit group only, then verify what remains unstaged.
5. If additional semantic groups remain after staging the first one, report them separately and leave them unstaged unless the user explicitly asks you to continue.

Repository-specific rules:

- Do not run `git commit`.
- Stage files with `git add` only.
- Never stage `claude-projects/`.
- Commit message proposals must not contain backticks.

Commit message rules:

- Keep commit messages concise: 1 or 2 sentences maximum.
- Be strictly factual. Do not exaggerate neatness, performance, maintainability, correctness, or scope.
- Describe what changed and, only when justified by the diff, why.
- Do not claim improvements that are not clearly supported by the changes.
- Do not rewrite, revert, or discard user changes unless explicitly asked.
- Do not combine unrelated edits just to reduce commit count.
- Prefer a small number of coherent commits over many tiny commits, but split distinct concerns when the diff supports it.

Before staging, show the planned commit list with a short title per group. Then stage the next group.

At the end, report:

- The group you staged
- The proposed commit message for that staged group
- Any changes intentionally left unstaged
- Any uncertainty you encountered while grouping changes

If the user supplied extra instructions as an argument, apply them as long as they do not conflict with the rules above.