# Branching and Releases

**Status: proposed.** Nothing described here exists yet: no GitHub App, no
rulesets, no `branch-policy` check, no Stage, Release, Manage staged release or
Resolve syncs actions. [Migration](#migration) lists what has to be built.

## Summary

There are three kinds of **mainline** branch:

- `dev` holds work toward the next release.
- `test/M.m.0` holds a **staged release**: a version being tested before it
  ships.
- `main` holds the **releases**: exactly one commit per published version.

People work on **side branches** and squash-merge them into their own
mainline. Moving changes from one mainline to another is never done by hand.
Four actions, run from the repository's Actions tab, make every such move,
every merge commit, and every tag:

| Action | Input | What it does |
|---|---|---|
| **Stage** | version, optional | moves `dev` into testing: stages a new release, or adds `dev`'s latest changes to the staged one |
| **Release** | version, optional | moves the staged release, or a hotfix, into `main`; tags it and publishes a GitHub Release |
| **Manage staged release** | operation, and a version | changes the staged release's version, or abandons it |
| **Resolve syncs** | none | lands sync branches whose conflicts a person has resolved |

![Side branches are squash-merged into their mainline; Stage moves dev into the staged release, and Release moves it or a hotfix into main; the bot syncs each move back.](images/Branching/branch-flow.svg)

How to read the diagram:

- **Rounded boxes** are mainlines; **dashed boxes** are side branches.
- **"squash"** arrows are people merging their own PRs.
- **"Stage"** and **"Release"** arrows are actions you run.
- **Dotted arrows** are syncs, which the bot makes after each move. After a
  hotfix, the sync goes into `dev` instead if nothing is staged.

Every change is checked before it lands, and no branch is ever force-pushed.
So mainlines are always in a valid state, and nobody's clone is ever rewritten.

## Terms

| Term | Meaning |
|---|---|
| **Release** | A version on `main`: one squash commit, its `vM.m.p` tag, and a GitHub Release. Nothing else is called a release. |
| **Staged release** | The open test branch `test/M.m.0`: a version being tested. At most one exists at a time. |
| **Release PR** | The PR from a staged release or a hotfix into `main`. Only the Release action merges it. |
| **Candidate** | A `vM.m.p-rcN` tag on a tested commit of a staged release or a hotfix branch. |
| **Release notes** | `docs/changelog/M.m.p.md`: one file per version, drafted by Stage, curated on the staged release, and published by Release. Together, the files are the changelog. |
| **Sync** | A merge commit, made by the bot, that carries a move between mainlines back to the branch it didn't land on. |
| **Mainline** | `main`, `dev`, or a test branch. |
| **Side branch** | Any other branch: `feature/`, `chore/`, `bugfix/`, `testfix/`, `hotfix/`, `sync/`. |

## Branches

| Branch | Started from | Merged into | By | Lifetime |
|---|---|---|---|---|
| `main` | — | — | — | permanent; the default branch |
| `dev` | — | the staged release | Stage | permanent |
| `test/M.m.0` | `dev`, by Stage | `main` | Release | staged while its release PR is open; frozen once released, deleted if abandoned |
| `feature/<name>`, `bugfix/<name>` | `dev` | `dev` | a person, squash | until merged |
| `chore/<name>` | `dev` or the staged release | the branch it started from | a person, squash | until merged |
| `testfix/<name>` | the staged release | the staged release | a person, squash | until merged |
| `hotfix/M.m.p` | `main` | `main` | Release | deleted once released |
| `sync/dev`, `sync/test-M.m.0` | the destination, by the bot | the destination | the bot, or Resolve syncs | until landed |

Rules that apply across the table:

- **At most one release is staged at a time.** A test branch is *staged* while
  its release PR `test/M.m.0 → main` is open. That single fact decides
  everything else:
  - **Released** (the PR was merged): the branch is *frozen*. It stays in
    place, because its sync into `dev` reads from it, but accepts no further
    changes.
  - **Abandoned** (the PR was closed without merging, by Manage staged
    release): the branch is deleted.
- **A side branch must start from the mainline it merges into.** Git doesn't
  record where a branch started, so this is checked at the PR: the PR may only
  bring in its own commits, and none that are already on another mainline. For
  example, a `feature/` started from `test/1.4.0` would carry that branch's
  test fixes into `dev`, so the PR is rejected with "rebase it onto dev".
  - **Commits already on the PR's own base never count**, so syncs that have
    already landed cause no false rejections. A branch started from `dev`
    after a sync, or updated by merging `dev` into it, still brings in only its
    own commits.
  - **Merging `main` into a side branch is rejected** while `main` holds a
    hotfix that hasn't reached that branch's base yet. The hotfix arrives
    through its sync.
  - **Sync PRs are exempt**, because carrying another mainline's commits is
    their purpose. Without the exemption, every sync would fail this check,
    and could never land if landing waited for its checks to pass. What a sync
    brings in is checked when it lands instead. Each of its merges must have
    the expected second parent (see [What a sync merges](#what-a-sync-merges)),
    and every commit it brings in is reachable from those parents. So that
    check already confines a sync to the move it records.

  These cases were checked in a throwaway repository:

  | PR | Result |
  |---|---|
  | `feature/` started from `dev` after a sync landed | passes |
  | `feature/` updated by merging `dev` | passes |
  | `testfix/` started from the staged release after a hotfix sync | passes |
  | `feature/` started from the staged release | rejected |
  | `feature/` that merged `main` while a hotfix hadn't reached `dev` | rejected |
  | A sync PR, `main → test` or `test → dev` | rejected by the rule, hence the exemption |

### What the history looks like

- **`main`** is linear: one commit per release, each tagged. `git log main` is
  the release list.
- **`dev` and test branches**, read with `git log --first-parent`, show one
  commit per squashed PR, per Stage, and per sync. Test branches also have the
  bot's release-notes commits. The syncs are the only merge commits, and the
  bot makes all of them.

## Stage

Stage moves `dev` into testing. Its one input, the version, says which release
to stage into. Leave it blank to mean the one already staged.

| You enter | A release is staged | Nothing is staged |
|---|---|---|
| blank | adds `dev` to it | error: "Nothing is staged. Enter a version to stage a new release." |
| the staged release's version | adds `dev` to it | — |
| another version | error: "1.4.0 is already staged. Release it, or change its version or abandon it with Manage staged release." | stages that version |

From the command line: `gh workflow run stage` or
`gh workflow run stage -f version=1.5.0`.

**Staging a new release** (`1.5.0`, nothing staged):

1. Checks:
   - the version has the form `M.m.0`;
   - it is higher than `main`'s version;
   - the tag `v1.5.0` does not already exist;
   - no sync into `dev` is open;
   - `main` is contained in `dev`, so every hotfix has reached it. This only
     fails after a staged release that had received hotfixes was abandoned. In
     that case Stage opens a sync from `main` into `dev` and stops; run Stage
     again once the sync has landed.
2. Creates `test/1.5.0` at `dev`'s tip.
3. Commits `version = "1.5.0"` in `pyproject.toml` on the new branch. The same
   commit deletes any leftover release-notes draft from an abandoned release.
4. Drafts the release notes, `docs/changelog/1.5.0.md`, in a second commit
   (see [Release notes](#release-notes)). Both commits are pushed together, so
   CI runs once and the first candidate is `-rc1`.
5. Opens the release PR `test/1.5.0 → main`, titled "Release 1.5.0". CI runs on
   it for the whole test period.
6. Syncs the version commit back into `dev`. The sync is built from the version
   commit alone, so the draft stays on the staged release.

**Adding to the staged release** (blank, or its version):

1. Checks: no sync into `dev` or into the staged release is open.
2. If `dev` doesn't yet contain the staged release's latest fixes, and
   squashing `dev` into it would conflict, Stage stops. It opens a sync from the
   staged release into `dev`, for a person to resolve. Run Stage again once
   that sync has landed; the squash then cannot conflict.
3. Squashes `dev` into the staged release, as "Stage dev into 1.4.0". The
   commit message records which `dev` commit was staged (a `Staged-from:`
   line).
4. Syncs the staged release back into `dev`. This sync carries everything on
   the staged release as it was before the squash, including the release-notes
   draft (see [Release notes](#release-notes) for why that is harmless).
5. Appends a block for the newly staged PRs to the end of the release notes, in
   a commit after the sync.

## Release

Release moves the staged release, or a hotfix, into `main`. Its input follows
Stage's pattern: blank means the staged release, and `1.4.1` means
`hotfix/1.4.1`. From the command line: `gh workflow run release`, or
`gh workflow run release -f version=1.4.1`.

1. Finds the release PR from that branch into `main`, and checks:
   - its required checks passed on its current head, which makes the head the
     latest candidate;
   - the version in the branch's `pyproject.toml` matches the branch name;
   - the version is higher than `main`'s;
   - the tag `vV` does not exist yet;
   - for a staged release: no sync into it is open, and `main` is already
     contained in it (so every hotfix has reached it);
   - `docs/changelog/V.md` exists, and it is the branch's only changelog file
     for a version without a release tag.
2. Writes the PR description:
   - the contents of `docs/changelog/V.md`;
   - a `Closes #N` line for each issue the release fixes, collected from the
     `Fixes #N` lines of the PRs it contains;
   - the candidate being promoted, for example "Promoted from v1.4.0-rc3".
3. Squash-merges the PR into `main` through the API, with that description as
   the commit message. Merging a PR into the default branch closes the issues
   it lists.
4. Tags the commit `vV`, builds the wheel and sdist, and publishes the GitHub
   Release, with `docs/changelog/V.md` as its text and both files attached.
5. Syncs the change back:
   - **a staged release** syncs into `dev`; the test branch is now frozen;
   - **a hotfix** syncs into the staged release, or into `dev` if nothing is
     staged; the hotfix branch is then deleted. A hotfix synced into a staged
     release reaches `dev` when that release ships.

## Manage staged release

Manage staged release changes the staged release's version, or abandons it.
Hotfixes are not managed this way: a hotfix branch belongs to the person making
it (see [Edge cases](#edge-cases) for renaming one). It takes two inputs:

| Input | Value |
|---|---|
| Operation | *Change version* (the default), or *Abandon* |
| Version | for *Change version*, the new version; for *Abandon*, the staged release's version, typed to confirm |

From the command line:
`gh workflow run manage-staged-release -f operation=change-version -f version=2.0.0`,
or `gh workflow run manage-staged-release -f operation=abandon -f version=1.4.0`.

### Change version

Use this when a staged release turns out to need a different version, for
example `1.4.0` becoming `2.0.0` because it breaks the API.

1. Checks:
   - the new version has the form `M.m.0`;
   - it differs from the current one and is higher than `main`'s;
   - its tag does not exist yet;
   - no sync into `dev` or into the staged release is open. Their instructions
     name the branch, which is about to change.
2. Creates `test/2.0.0` at the staged release's tip. It commits
   `version = "2.0.0"` on it, and in the same commit renames the release notes
   to `docs/changelog/2.0.0.md`, keeping their content.
3. Moves every open PR into the staged release over to `test/2.0.0`. Each
   one's version check now asks for `2.0.0`, so it must be updated from the new
   branch before it can merge.
4. Closes the release PR "Release 1.4.0", and opens "Release 2.0.0".
5. Deletes `test/1.4.0`.

The old version's candidate tags are kept, and the new version's candidates
start again from `-rc1`. `dev` picks up the new version when the release syncs
back.

### Abandon

Use this when a staged release will not ship. If its changes should ship
later, run Stage again from `dev`.

1. Checks that the version names the staged release.
2. Closes, each with a comment saying the staged release was abandoned:
   - the release PR;
   - every open PR into the staged release;
   - every open sync into or out of it. That is its own sync branch, plus
     `dev`'s, which carries only merges from the staged release while one
     exists.
3. Deletes the test branch.

**Nothing is synced.** Abandon discards everything that exists only on the test
branch: the `testfix/` and `chore/` PRs squashed into it, and whatever its
unlanded syncs held. The run summary, and the comment on the release PR, list
those PRs, so any of them can be redone on `dev`. Work that came from `dev` is
unaffected, since `dev` still has it.

Hotfixes synced into the abandoned staged release are not lost. They are on
`main`, and the next Stage brings them into `dev` before it stages a release.

The release-notes draft goes with the branch. If Stage had added `dev` during
testing, `dev` also has a copy, which the next new release's version commit
deletes (see [Release notes](#release-notes)). Either way, the next draft
starts clean and lists every unreleased PR. Curation of the abandoned draft is
lost.

Abandon is never blocked. Candidate tags of the abandoned version are kept. If
the same version is staged again, its candidates continue from the next number.

## Syncs

A **sync** carries a move back to the branch it didn't land on. Each sync is a
real merge commit, built by the bot on a `sync/` branch:

1. The bot opens the sync as a PR.
2. If the merge was clean and CI passes, the bot lands it on its own.
3. Otherwise the PR stays open with the `open-sync` label:
   - a person resolves the conflict on the sync branch, using ordinary git;
     the PR's instructions give the exact commands;
   - then they run **Resolve syncs**, which checks the result and lands it
     (see [Safety checks](#safety-checks)).

People never merge sync PRs themselves.

### What a sync merges

A squash commit shares no history with the branch it came from. If a sync
merged the squash commit itself, git would compare against an old common
ancestor. Every line edited again after the Stage would then show as a
conflict, and these conflicts grow with every cycle. So each sync takes its
**content** from a branch that still shares history, and records the **squash
commit** as its second parent, so the next move starts from it:

| Move | Sync goes into | Content comes from | Second parent |
|---|---|---|---|
| Stage a new release | `dev` | the version commit on the new test branch | the same |
| Stage, adding to the staged release | `dev` | the staged release as it was before the squash | the Stage squash commit |
| Release the staged release | `dev` | the frozen test branch | the release commit on `main` |
| Release a hotfix | the staged release, else `dev` | `main` | the same |
| Stage, after an abandoned staged release that had received hotfixes | `dev` | `main` | the same |

These were checked in a throwaway repository, over three release cycles in
which `dev` kept re-editing lines that had already been staged:

| Variant | Result |
|---|---|
| Syncing the squash commit itself | the first release sync conflicted; by the second cycle the release squashes conflicted too |
| As in the table | every Stage, release and sync merged cleanly; `main` stayed linear |

When conflicts arise, the instructions are plain git, for example
`git merge origin/test/1.4.0`. Resolve syncs adds the second parent when it
lands the sync.

### When an open sync blocks an action

An action is blocked only by open syncs **into** a branch it reads from or
writes to:

| Action | Reads → writes | Blocked by open syncs into |
|---|---|---|
| Stage, new release | `dev` → new test branch | `dev` |
| Stage, adding to the staged release | `dev` → the staged release | `dev` and the staged release |
| Release, staged release | the staged release → `main` | the staged release |
| Release, hotfix | `hotfix/*` → `main` | nothing: no sync ever goes into `main` |
| Manage staged release, change version | the staged release → new test branch | `dev` and the staged release |
| Manage staged release, abandon | deletes the staged release | nothing: it closes the syncs it would wait on |

So a hotfix can always be released. A staged release waits only for hotfixes
into it. A sync into `dev` that conflicts holds up Stage, but never a release.

### One open sync per destination

Sometimes a second sync heads for a branch that already has one open. For
example, a hotfix ships while the release sync into `dev` is still unresolved.
The bot then appends the new merge to the open sync branch, after the existing
one. It does not open a second PR, which would merge from an out-of-date
starting point and repeat the first PR's conflicts. The `open-sync` PR lists
its merges as a checklist. If an earlier merge still has unresolved conflicts,
Resolve syncs adds the later ones once it is resolved.

Appending was checked in the same throwaway repository, on the awkward case:
`dev`'s sync branch received, in turn, a Stage sync, then a release sync, then
a hotfix sync. All three merged cleanly, the content was correct, and the next
release cycle was clean.

### Why this cannot deadlock

Actions wait only on syncs, and syncs wait only on a person resolving them.
Resolve syncs is never blocked, and no sync needs any action to run first.
Nothing ever waits on an action, so the waiting can never form a cycle.

### Landing a sync

A sync branch is built on the destination's tip:

- **If the destination hasn't moved** when the sync lands, landing
  fast-forwards it.
- **If people have squashed side branches into it meanwhile** (`dev` moves
  often), the bot merges the sync branch into the current tip. That adds one
  merge commit.

In both cases `git log --first-parent` shows the sync as a single commit, and
nothing is force-pushed.

## Versioning

**The version is the static `version` in `pyproject.toml`.**

- `main` carries its latest release.
- Each staged release `test/V` carries `V`. Stage sets it when it stages the
  release, and Manage staged release when it changes the version.
- `dev` follows automatically: syncs bring in each newly staged version and
  each hotfix version.
- A hotfix's version is set by the person making it.

**Form of the version:**

- Staged releases are always `M.m.0`.
- Patch versions come only from hotfixes, named `hotfix/M.m.p` with `p ≥ 1`.
- Choosing between a minor and a major release happens when you type the
  version into Stage. It can be changed later with Manage staged release.

**Checks.** The failure messages tell the person exactly what to set; nothing
increments versions automatically. Every check, with what it prevents, is in
[Safety checks](#safety-checks).

| Where | Check |
|---|---|
| Branch names | `test/M.m.0`; `hotfix/M.m.p` with `p ≥ 1` |
| PR into `test/V` | the incoming branch's `pyproject.toml` says `V` |
| Stage, Release, Manage staged release (change version) | the version is higher than `main`'s, and its tag does not exist yet |
| Release | the branch's `pyproject.toml` matches its name |

**The one version conflict.** A hotfix synced into a staged release always
conflicts on the version line, because both sides changed it. The staged
release's version is the right answer, and the check on PRs into `test/V` says
so. When that line is the only conflict, the bot resolves it itself.

**Semantic versioning, with 0.x semantics until the API settles.**

- In `0.x`, a minor release may break the API; a patch release never does.
- Breaking changes are labelled `breaking`, and the release-notes draft lists
  them first.

The versions are plain PEP 440 versions (`1.4.0`, `1.4.1`). The `v` on tags is
a git convention that packaging tools strip. Candidate tags such as
`v1.4.0-rc2` follow PEP 440 too: tools read them as `1.4.0rc2`, which sorts
before `1.4.0`. Candidates are tags only; `pyproject.toml` is not changed for
them (see [Release candidates](#release-candidates)). There is no `.devN`
suffix: nothing is published from `dev`, and `dev`'s version always names a
version that has been released or staged.

## Releases

A release is a squash commit on `main`, the `vM.m.p` tag on it, and a GitHub
Release.

- **Release notes** live in `docs/changelog/`; see
  [Release notes](#release-notes).
- **Squash commit messages** are the PR's title and description (a repository
  setting). So each `dev` commit carries its PR's `Fixes #N` lines, which
  Release turns into `Closes #N`.
- **Artifacts.** Every release PR has a required check that builds the sdist
  and wheel, installs the wheel into a clean environment, and imports `PyPR`.
  Release builds them again and attaches them to the GitHub Release.
- **Installing a release:**
  `pip install "git+https://github.com/DavidCGordon/Product-Registers@v0.1.0"`.
  A candidate can be installed the same way, with `@v0.1.0-rc2`.
- **No PyPI for now.** Releases are installed from git tags. If that ever
  changes, note that the name `PyPR` on PyPI belongs to an unrelated project,
  so publishing would need a different distribution name (the import name can
  stay `PyPR`), and Release would publish with PyPI trusted publishing, which
  needs no stored token.
- **Documentation** is versioned along with the code; see
  [Documentation versions](#documentation-versions).

### Release notes

Each release's notes are a file, `docs/changelog/M.m.p.md`. Together the files
are the changelog. The docs' Changelog page lists them newest first, and each
GitHub Release shows its own version's file.

**One file per version.** Only one branch ever writes a given version's file:
the staged release, or the hotfix branch. So the files never conflict. A single
`CHANGELOG.md` would conflict constantly, because every release would edit the
same lines at the top.

**What a draft contains.** The template is a heading and an empty Highlights
section for you to write. Each Stage then appends one block to the end of the
file. A block has one line per PR title for the `feature/` and `bugfix/` PRs
that Stage brought in, grouped under Breaking changes (PRs labelled
`breaking`), Features, and Fixes. `chore/` and `testfix/` PRs are left out,
since a test fix usually repairs something that was never released; add any
that matter by hand.

**How each action touches the draft:**

| Action | What it does to `docs/changelog/` |
|---|---|
| Stage a new release | deletes any changelog file for a version without a release tag (a leftover from an abandoned release), then writes the template for the new version, with a block of every PR since the previous release |
| Stage, adding to the staged release | appends a block of the newly staged PRs to the end of the draft |
| Change version | renames the draft to the new version, and updates its heading |
| Abandon | nothing; the draft goes with the test branch |

**The draft's content is never read.** Stage only writes a fresh file or
appends to the end of one, so your edits are never touched, and the file can
take any shape. It knows which PRs are new from git, not from the file: each
Stage's squash records the `dev` commit it staged, in its `Staged-from:` line.
The new PRs are those between the previous Stage's `Staged-from:` commit and
this one's, or since the previous release when staging a new release.

**Leftovers from an abandoned release.** If Stage had added `dev` before the
abandon, `dev` holds a copy of the abandoned draft. The next new release's
version commit deletes it, and that commit's sync carries the deletion into
`dev`. The new draft then starts clean, listing every unreleased PR once,
including those of the abandoned release. Curation from the abandoned release
is lost.

Deleting the copy this way is safe. A deletion causes trouble only when a
later merge involves a branch that still has the old file, as explained under
"Where the draft lives" below. After an abandon, no such branch is
left: the abandoned test branch and its syncs are gone, and the new test
branch starts after the deletion.

**Curating.** Edit the draft in a `testfix/` or `chore/` PR into the staged
release: write the highlights, merge the blocks, reword lines, drop noise.
Release publishes the file exactly as it stands, so an uncurated block ships
as written.

**Where the draft lives.** Stage commits the draft after the version commit, and
syncs only the version commit. So unless more of `dev` is staged, `dev` first
sees the notes in the release sync, final.

When Stage adds `dev` to a staged release, its sync carries the staged
release's content into `dev`, draft included. That cannot be avoided. A sync
that left the file out would be recorded as `dev` deleting it, and later
merges would delete it from the staged release, or conflict. So a draft can sit
on `dev`. It is harmless there:

- `branch-policy` rejects `dev` PRs that edit the changelog file of any
  version without a release tag, so only the staged release changes the draft,
  and later syncs simply update the copy.
- `dev`'s docs are never hosted, so no user sees it.

**Hotfixes.** The author writes `docs/changelog/M.m.p.md` on the hotfix branch.
`branch-policy` requires it on hotfix PRs into `main`, and its failure message
includes the template. If the hotfix's version changes, the author renames the
file to match.

**Publishing.** Release checks that the file exists, and that no other
unreleased version's file is present. It then uses the file as the GitHub
Release's text and in the release commit's message. So `main`, and every
hosted docs build, contains only final notes.

**The Changelog page** is generated by `docs/conf.py` from `docs/changelog/`,
sorted by version number, so `0.10.0` comes after `0.9.0`. Each version's docs
show the history up to that version, and candidate docs show the draft as it
currently stands.

### Documentation versions

The docs site is not hosted yet; `.readthedocs.yaml` is in place for Read the
Docs. Hosted docs are built only from version tags. Unreleased work, on `dev`
or anywhere else, is never published.

**What users see:**

| Hosted version | Built from | Shown |
|---|---|---|
| `stable` | the newest release tag | the site's front page |
| `v0.1.0`, `v0.2.0`, … | each release tag | listed in the version menu |
| `v0.2.0-rc1`, `-rc2`, … (optional) | each candidate tag, on a staged release or a hotfix branch | hidden: reachable only by link, for reviewing a release's docs before it ships |

So users see the docs of the newest release by default, and can switch to any
earlier release. The docs change only when a release ships.

**How docs move through the flow:**

1. Docs change on side branches like any other change, and land on `dev`. To
   see them rendered, build locally (see the README).
2. Stage puts them on the staged release. Each candidate tag builds a hidden
   version, and the bot links it from the release PR.
3. Release tags `vM.m.0`. Read the Docs builds it, lists it, and moves
   `stable` to it.
4. A hotfix works the same way: candidate docs under `vM.m.p-rcN`, and `stable`
   moves when Release tags `vM.m.p`.

**What each build shows:**

- **The title names the tag:** "PyPR 1.4.0" for a release, "PyPR 1.4.0-rc2"
  for a candidate. A build of `main` is titled with its version too, since
  `main` is always exactly the newest release. Every other build, local builds
  included, is titled "PyPR (unreleased)". No branch name starts with `v`, so
  `docs/conf.py` can tell a tag from a branch by its first letter.
- **Links point at the commit being built.** Links that leave `docs/`, to
  tests, the README, or source files, point at the exact commit the build was
  made from. A release's docs therefore never link to newer code, or to
  `master`.
- **The API reference** is generated from the docstrings on every build, so it
  always matches its version.
- **The Changelog page** shows the release notes up to that version (see
  [Release notes](#release-notes)).

**Builds are automatic.** Read the Docs is connected to the repository once,
and from then on hears about every pushed tag, including tags the bot pushes.
The automation rules below decide which tags to build. Nothing changes per
release, and no file in the repository lists the versions.

**A release's docs are fixed once it ships.** Each version is built from its
tag, and tags cannot be moved (see the rulesets), so a published version's docs
cannot be edited in place. That guarantees the docs for `v1.4.0` describe
exactly the code in `v1.4.0`. To fix a mistake in published docs:

- **Usually:** fix it on `dev`. It ships with the next release, when `stable`
  moves to it.
- **When it can't wait** (wrong instructions, a misleading API description): a
  hotfix that changes only docs, for example `hotfix/1.4.1`. It goes through
  candidates and Release like any hotfix, and `stable` moves to `v1.4.1` as soon
  as it ships. A patch release for a documentation fix is normal under semantic
  versioning.

Moving a tag, or keeping a separate branch for docs fixes, would make docs
editable after release. But it would break the link between a version's docs
and its code, and add a mainline outside the release flow, so this design does
neither.

**Read the Docs settings,** set once when the project is created:

| Setting | Value | Why |
|---|---|---|
| Default version, where the site's root goes | `stable` | users land on the newest release |
| `stable` | automatic | Read the Docs points it at the highest release tag, and skips pre-releases such as `-rcN` |
| `latest` | hidden from the version menu | Read the Docs always keeps a `latest` version, built from GitHub's default branch. That is `main`, the newest release, so hiding it loses nothing and keeps `dev` out entirely |
| Automation rule | activate tags matching `^v\d+\.\d+\.\d+$` | every release is built and listed. The predefined SemVer rule would also match the `-rcN` tags, so a custom pattern is used |
| Automation rules (optional, for candidate docs) | for tags matching `^v\d+\.\d+\.\d+-rc\d+$`: activate, and hide | candidates are built but never listed |
| Pull request builds | off | they would publish previews of unreleased work |

**A release's docs must build.** `tests.yml` runs `sphinx-build -W` on every
PR (warnings count as errors; the build has none today), and Release requires
it to have passed on the release PR. A release tag can then never fail to build
on Read the Docs.

### Release candidates

A candidate is a commit that passed CI and may become a release. Whenever the
tip of a staged release or of a hotfix branch changes, and CI passes on the new
tip, the bot tags it with that branch's version and the next candidate number:
`v1.4.0-rcN` on `test/1.4.0`, `v1.4.1-rcN` on `hotfix/1.4.1`. A tip that fails
CI gets no tag.

The tip of a staged release changes whenever:

- Stage stages it, with its version commit;
- a `testfix/` or `chore/` PR is squashed in;
- Stage adds `dev` to it;
- a sync lands;
- Manage staged release changes its version.

A hotfix branch's tip changes whenever its author pushes to it.

`N` is one more than the highest existing candidate for that version, rather
than a count of commits. That keeps numbers unique when one push carries
several commits, and when an abandoned version is staged again. A version
change starts again from `-rc1` under the new version. This applies equally to
a staged release changed with Manage staged release, and to a hotfix branch
renamed by hand.

- **Releasing promotes the latest candidate.** Release requires the branch's
  tip to have passed CI, so that tip is always a tagged candidate, and the
  release notes name it.
- **Candidates of a hotfix outlive its branch.** The hotfix branch is deleted
  once released, but its candidate tags remain, and keep their commits
  reachable.
- **Candidates are tags, not versions.** The package installed from
  `@v1.4.0-rc2` still reports version `1.4.0`, because `pyproject.toml` is not
  changed for candidates. Rewriting it on every commit would add a bot commit
  to each change, and would break the check on PRs into `test/V`. The tag
  identifies the exact build being tested, which is what candidates are for.

## Enforcement and automation

The repository is public, so rulesets are available on the free plan.

### Repository settings

| Setting | Value |
|---|---|
| Default branch | `main` |
| Merge buttons | squash only: merge commits and rebase merging off, since the bot lands its merge commits by pushing |
| Squash commit message | pull request title and description |
| Automatically delete head branches | on |

### Rulesets

| Ruleset | Targets | Rules | Bypass |
|---|---|---|---|
| Main | `main` | *Restrict updates* (only bypass actors may change it), no deletion or force push, linear history | the App |
| Dev and test | `dev`, `test/*` | PR required (0 approvals), required checks `tests-passed`, `branch-policy`, `manual-merge`; no deletion or force push | the App |
| Bot-made branches | `test/*`, `sync/*` | restrict creations | the App |
| Branch names | every branch except `main`, `dev`, and the patterns `test/*`, `sync/*`, `feature/*`, `chore/*`, `bugfix/*`, `testfix/*`, `hotfix/*` | restrict creations | the App |
| Version tags | `v*`, which includes the `-rcN` candidates | restrict creation, update and deletion | the App |

A release PR into `main` shows its checks passing and its merge box blocked:
only the App can update `main`. The bot comments on it with how to release it.

### Checks

- **`branch-policy`** (required). It runs on `pull_request_target`, which
  runs the workflow as it exists on the base branch, so a PR cannot edit the
  check that judges it. It reads git refs, and parses `pyproject.toml` from the
  PR as data; it never runs code from the PR. It checks:
  - **allowed sources:**

    | Into | Accepts |
    |---|---|
    | `dev` | `feature/`, `chore/`, `bugfix/` |
    | `test/V` | `testfix/`, `chore/` |
    | `main` | `test/V`, `hotfix/V` |
    | any | the matching `sync/` branch |

  - the branch-start rule, for people's side branches; sync PRs are exempt;
  - the version checks;
  - that a test branch accepts PRs only while it is staged;
  - that a PR into `dev` edits no changelog file of a version without a
    release tag, and that a hotfix PR into `main` includes its own
    `docs/changelog/M.m.p.md`.
- **`manual-merge`** (required). The App reports it:
  - **passed** for side-branch PRs;
  - **"action required: landed by Resolve syncs"** for sync PRs.

  This keeps the squash button disabled on sync PRs. It has to be reported
  explicitly: GitHub counts a required check whose job was *skipped* as
  passed.
- **`tests-passed`**, from `tests.yml`. Its triggers change to pushes to
  `main`, `dev`, `test/**` and `hotfix/**`, and PRs into any mainline. The
  hotfix pushes are what give hotfix candidates a tested tip to tag.
- **`docs-build`**, a new job in `tests.yml`: `sphinx-build -W` on every PR.
  Release requires it on the release PR.

### Safety checks

Every check in the system, where it runs, and what it prevents. Each failure
stops the action or the merge with a message saying what to do. None of these
checks change anything themselves.

| Check | Where | Prevents |
|---|---|---|
| The source branch's prefix is allowed for the base | `branch-policy` | work landing on the wrong mainline |
| A PR from a person's side branch brings in only its own commits; sync PRs are exempt | `branch-policy` | a side branch started from the wrong mainline carrying that mainline's commits along, for example a staged release's fixes into `dev` |
| A test branch accepts PRs only while it is staged | `branch-policy` | changes to a frozen or abandoned test branch |
| Sync PRs cannot be merged by people | `manual-merge` | a squashed sync, which would lose its link to the history it merges |
| A PR cannot edit the check that judges it | `branch-policy` runs on `pull_request_target`, and never runs PR code | a PR weakening its own check |
| Branch names are `test/M.m.0` and `hotfix/M.m.p` with `p ≥ 1` | `branch-policy`, Stage, Manage staged release | malformed or out-of-place versions |
| A PR into `test/V` carries version `V` | `branch-policy` | the staged release's version drifting; also tells whoever resolves a hotfix sync which version line to keep |
| The version is higher than `main`'s | Stage, Release, Manage staged release (change version) | `main`'s version going backwards, for example a `1.4.1` hotfix shipped after `1.5.0` |
| The version's tag does not exist yet | Stage, Release, Manage staged release (change version) | reusing a version that was already released |
| The branch's `pyproject.toml` matches its name | Release | a release whose package version differs from its tag |
| The release PR's checks passed on its current head | Release | shipping a commit that is not the latest tested candidate |
| The docs build with no warnings | `docs-build`, required by Release | a release whose docs fail to build on Read the Docs |
| A PR into `dev` edits no changelog file of an unreleased version | `branch-policy` | the copy of a draft on `dev` diverging from the staged release's, which would conflict at the next sync |
| A hotfix PR into `main` includes `docs/changelog/M.m.p.md` for its version | `branch-policy` | a hotfix release with no notes |
| The release has its own changelog file, and no other unreleased version's | Release | a release without notes, or a stale draft reaching `main` and the hosted docs |
| `main` is contained in the staged release | Release (staged release) | shipping a release that lacks a hotfix |
| `main` is contained in `dev` | Stage (new release) | staging a release that lacks a hotfix, after an abandoned staged release had received one |
| No open sync into a branch the action reads or writes | Stage, Release, Manage staged release (change version) | building on a branch that is missing changes, and spurious conflicts later |
| Squashing `dev` into the staged release must not conflict | Stage (adding) | resolving conflicts on a mainline; Stage opens a sync for the conflict instead |
| The version is typed to confirm | Manage staged release (abandon) | abandoning a staged release by accident |
| The sync branch has no conflict markers, its merges have the expected parents, and CI passed on its head | Resolve syncs, and automatic landing | landing a half-resolved or wrongly built sync, or one that carries anything beyond the move it records |
| The version line is resolved automatically only when it is the only conflict | the sync builder | the bot deciding a real conflict |
| A candidate is tagged only after CI passes on it | `candidates.yml` | an untested candidate |
| Actions run one at a time, never force-push, and recompute if a branch moved | all actions | races between actions and people, and rewritten history |

### The GitHub App

The bot is a GitHub App: an identity registered under your account settings,
installed on this repository, with read/write access to contents and pull
requests. Its ID and private key are stored as repository secrets, and the
workflows exchange them for a short-lived token with
`actions/create-github-app-token`. Its actions show as `<name>[bot]`.

The default `GITHUB_TOKEN` cannot do this job, for two reasons:

- **Nothing it does triggers other workflows.** GitHub does this to prevent
  loops. Sync PRs opened with it would never run CI, and their required checks
  would never report.
- **It cannot be put on a ruleset's bypass list**, and the App must bypass
  rulesets to update `main`, create test and sync branches and tags, and land
  merge commits.

### Workflows

| Workflow | Trigger | Purpose |
|---|---|---|
| `branch-policy.yml` | PRs into a mainline | the two checks above, and the "ready" comments |
| `stage.yml` | manual; input: version | Stage |
| `release.yml` | manual; input: version | Release |
| `manage-staged-release.yml` | manual; inputs: operation, version | Manage staged release |
| `resolve-syncs.yml` | manual; also when `Tests` completes on a `sync/` branch | Resolve syncs; lands bot-made clean syncs automatically once CI passes |
| `candidates.yml` | when `Tests` completes on a `test/` or `hotfix/` branch | tags a passing tip as the next `-rcN` |
| `tests.yml` | push and PR | unchanged apart from its branch list |

The four actions share one `concurrency` group, so only one of them runs at a
time. The sync-building logic lives in one script under `.github/scripts/`,
shared by all four. The actions never force-push. If a person's squash moves a
branch while an action runs, the push is refused, and the action recomputes on
the new tip.

Manually run workflows use the workflow file from `main` by default, so the
release machinery that runs is always the released version. A change to it on
`dev` takes effect at the next release.

## Day to day

- **A feature or fix:**
  1. Branch from `dev` as `feature/<name>` or `bugfix/<name>`.
  2. Open a PR into `dev`. GitHub suggests `main`, the default branch, so
     change the base.
  3. Squash-merge it once its checks pass.
- **Staging a release:** run Stage with the version, for example `1.5.0`. The
  release PR opens, and CI runs on it.
- **Fixing something in the staged release:**
  1. Branch from `test/1.5.0` as `testfix/<name>` (or `chore/<name>`).
  2. Open a PR into the staged release, and squash-merge it.

  The fix reaches `dev` when the release ships.
- **Adding newer `dev` work to the staged release:** run Stage with the field
  blank. It appends the newly staged PRs to the end of the release notes, for
  you to sort into place.
- **Curating the release notes:** edit `docs/changelog/1.5.0.md` in a
  `chore/` PR into the staged release. Write the highlights, reword lines, and
  drop noise.
- **Trying the current candidate:** install from the newest `-rcN` tag.
- **Shipping:** run Release with the field blank.
- **Changing the staged release's version:** run Manage staged release,
  choosing *Change version*, with the new version.
- **Dropping the staged release:** run Manage staged release, choosing
  *Abandon*, with its version. Its test fixes are discarded; the run summary
  lists them.
- **A hotfix:**
  1. Branch from `main` as `hotfix/1.5.1`, set `version = "1.5.1"`, and make
     the fix.
  2. Write its notes in `docs/changelog/1.5.1.md`.
  3. Push. Each push that passes CI becomes a candidate, `v1.5.1-rcN`.
  4. Open a PR into `main`.
  5. Run Release with `1.5.1`.
- **A sync conflict:**
  1. Check out the `open-sync` PR's branch and run the merge command it gives.
  2. Resolve, commit, and push.
  3. Run Resolve syncs.

## Edge cases

- **A hotfix ships while a release is staged.** Its sync goes into the staged
  release, and Release for the staged release waits until the sync has landed.
  This is automatic when the merge is clean, apart from the version line, which
  the bot resolves.
- **A release ships while a hotfix PR is open.** The hotfix's version may now
  be lower than `main`'s (it was `1.4.1`, and `main` is at `1.5.0`). Release
  refuses it. Rename the branch (GitHub keeps the PR attached), set the new
  version, and update the branch from `main`. Its candidates start again from
  `-rc1` under the new version.
- **Two hotfixes at once.** The branch name is the version, so the second one
  takes the next patch number after the first is released.
- **A fix made on `dev` that the staged release also needs.** Either run Stage
  to bring `dev` in, or make the fix as a `testfix/` instead.
- **Only the latest release line gets hotfixes.** There is one `main`, so
  `1.4.x` stops receiving fixes once `1.5.0` ships.
- **A staged release that turns out to be major.** Change its version with
  Manage staged release. Its test fixes carry over; its candidate numbering
  restarts under the new version.
- **Open PRs into an abandoned staged release** are closed. Their branches
  started from the test branch, so they carry its commits, which `dev` doesn't
  have. Re-opening the same branch against `dev` therefore fails the
  branch-start rule. To keep the work, start a new branch from `dev` and
  cherry-pick the PR's own commits onto it.
- **Dependabot and other bots** create branches with their own prefixes. Add
  their pattern to the allowed list, or leave them disabled.
- **Required reviews.** A sole maintainer cannot approve their own PR, so
  require zero approvals; the required checks still apply. Raise the count when
  there are other maintainers.

## Migration

`master` is left exactly as it is for now, so that clients pulling it are not
disrupted, and will be deleted later. `main` is created alongside it.

1. **Build the tooling** in one PR into `dev`:
   - the shared scripts and their tests;
   - the workflows and the rulesets, as JSON;
   - the Changelog page in the docs;
   - a PR template;
   - `__version__`, read with `importlib.metadata.version("PyPR")`, so
     `pyproject.toml` stays the only place the version is written.
2. **Rehearse** every flow in this document in a throwaway public repository,
   then delete it.
3. **Configure this repository:**
   - create and install the GitHub App;
   - create `main` at `master`'s tip, so it shares `master`'s history;
   - apply the repository settings and rulesets;
   - create the `open-sync` and `breaking` labels;
   - temporarily make `dev` the default branch. GitHub only runs a
     manually triggered workflow whose file is on the default branch, and
     `main` has none until the first release.
4. **First release:** run Stage with `0.1.0`, test the candidates, then run
   Release.
   - There is no previous release, so the draft starts from the template with
     no PR lines. Write the notes by hand in a `chore/` PR.
   - The release PR's `branch-policy` check is missing this once, because it
     runs from `main`'s copy of the workflow. Release's own checks stand in.
5. **Cut over:** make `main` the default branch, and stop running CI on
   `master`. Later, announce a date, then delete `master`.
