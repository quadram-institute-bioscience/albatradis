import argparse
import semantic_version


def main():
    parser = argparse.ArgumentParser(description=
                                     'Bumps the semantic version of this package.')
    parser.add_argument('--mode', type=str, required=True, default="patch",
                        help='Whether to bump (major, minor or patch) number (default: patch)')

    args = parser.parse_args()

    mode = args.mode
    if not mode:
        raise ValueError("User must specify which semantic version mode to use")

    version = open('VERSION').read().strip()
    oldsemver = semantic_version.Version(version)

    if mode == "major":
        newsemver = oldsemver.next_major()
    elif mode == "minor":
        newsemver = oldsemver.next_minor()
    elif mode == "patch":
        newsemver = oldsemver.next_patch()
    else:
        raise ValueError("Unrecognised semantic version mode")
    print("Updating Albatradis version from", oldsemver, "to", newsemver)
    open('VERSION', 'w').write(str(newsemver))


if __name__ == '__main__':
    main()