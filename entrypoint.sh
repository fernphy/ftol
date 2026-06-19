#!/bin/bash
set -e

HOST_UID=${HOST_UID:-0}
HOST_GID=${HOST_GID:-0}

# If no remapping requested, run as root (preserves current behaviour)
if [ "$HOST_UID" -eq 0 ]; then
    exec "$@"
fi

# Create group matching host GID if not already present
if ! getent group "$HOST_GID" >/dev/null 2>&1; then
    groupadd -g "$HOST_GID" hostgroup
fi

# Create user matching host UID if not already present
if ! getent passwd "$HOST_UID" >/dev/null 2>&1; then
    useradd -u "$HOST_UID" -g "$HOST_GID" -m -s /bin/bash analyst
fi

USER_NAME=$(getent passwd "$HOST_UID" | cut -d: -f1)

# Transfer /renv ownership so R package installs work without root
chown -R "$HOST_UID:$HOST_GID" /renv

exec gosu "$USER_NAME" "$@"
