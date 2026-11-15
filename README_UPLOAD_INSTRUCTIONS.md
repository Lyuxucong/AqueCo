# Upload the English AqueCo Manual to GitHub

This ZIP contains the `docs/` folder, plus this instruction file.

## Recommended integration

1. Back up the existing AqueCo repository.
2. Copy all files from this ZIP's `docs/` directory into the package repository's `docs/` directory.
3. **Merge**, do not delete, existing documents such as `QUICKSTART_CN.md`, `HELP_SYSTEM_CN.md`, or `AUTHOR_UPDATE_CN.md`.
4. Add this line to the package's root `README.md`:

   `📘 [Complete AqueCo User Manual (English)](docs/README.md)`

5. Commit and push the changes:

   ```bash
   git add docs README.md
   git commit -m "Add comprehensive English user manual"
   git push
   ```

No Quarto, GitHub Pages, or tutorial website is needed.

## Included material

- 23 English chapter pages and a `docs/README.md` table of contents.
- 19 figure placeholder PNGs in `docs/images/`.
- Examples and signatures grounded in AqueCo 0.1.0 source documentation.
- A 64-function exported API reference.

## Important

This is a source-grounded documentation draft, not a claim that all examples were executed. Verify each chapter in R and replace figure placeholders with actual results before describing the tutorial as fully validated.
