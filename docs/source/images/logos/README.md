# refineGEMs logo files

This directory contains the canonical refineGEMs logo assets. Choose a layout
according to the available space:

- `refineGEMs_full_*`: the complete, square logo. Use it when the logo can be
  shown at a size large enough for all details and text to remain legible.
- `refineGEMs_icon_*`: the symbol without the wordmark. Use it for compact,
  square placements such as application icons and avatars.
- `refineGEMs_icon_badge.png`: a 32×32 optimized copy of the transparent icon
  for embedding in badges. Use this asset instead of embedding a full-size
  canonical icon.
- `refineGEMs_text_*`: the horizontal wordmark. Use it for README headers,
  documentation headers, and other wide placements.

Each layout is available with three background variants:

- `*_nb.png`: no background (transparent)
- `*_wb.png`: white background
- `*_bb.png`: black background

Prefer the transparent variant when the surrounding background is controlled
and provides sufficient contrast. Use a background variant when the rendering
environment does not handle transparency reliably or when a consistent
background is required.

## Compatibility file

`../refineGEMs_logo.png` is an exact copy of `refineGEMs_full_nb.png`. It keeps
absolute image links in older refineGEMs release READMEs working after the
previous logo was retired. The compatibility filename does not identify a
separate logo and should not be used in new material; use one of the canonical
files in this directory instead.

Do not remove or replace the compatibility file with the retired artwork.

## Release badge

The GitHub release badge in the project README uses
`refineGEMs_icon_badge.png` as a custom Shields logo. Shields requires custom
logos to be embedded as a Base64 data URI. To refresh the badge after changing
the canonical icon:

```shell
sips --resampleHeightWidth 32 32 refineGEMs_icon_nb.png \
  --out refineGEMs_icon_badge.png
base64 < refineGEMs_icon_badge.png | tr -d '\n'
```

Insert the resulting value into this reusable Markdown pattern:

```md
![GitHub release (with filter)](https://img.shields.io/github/v/release/draeger-lab/refinegems?logo=data:image/png;base64,BASE64_DATA&label=refineGEMs&color=B4A069&style=flat-square)
```

For licensing and attribution, see [LOGO_LICENSE.md](LOGO_LICENSE.md). For
colours, sizing, spacing, and modification guidance, see
[the brand usage guide](../../brand.rst).
