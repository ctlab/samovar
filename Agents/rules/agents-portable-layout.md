# Portable SAMOVAR harness

**Conditional:** apply only when creating, changing or packaging SAMOVAR agent
skills, rules or hooks under `Agents/`.

Maintain IDE-neutral `SKILL.md` instructions, ordinary Markdown rules and small
executable helpers. Put IDE discovery/installation differences in the installer,
not in copies of each skill. No host paths, local task-manager clients, secrets,
model names or assumptions about a particular IDE's tools belong in the bundle.

Every rule must name its SAMOVAR-specific activation condition. Skills must link
the rules they need. Bundle these references during installation; users must
not need a second rules installation. Do not install unconditional global
`AGENTS.md`/`CLAUDE.md` policy for SAMOVAR usage.

Only skills listed in `bundle.json` are distributed. Existing `.agents/skills`
are checkout-maintenance instructions and are not distribution input. Add an
IDE hook only for a concrete requirement, with a tested adapter and explicit
installation consent; do not write lifecycle settings as an incidental install.
