#!/bin/bash
set -euxo pipefail

dnf update -y
dnf install -y git docker perl
systemctl enable --now docker

if [ ! -f /swapfile ]; then
  fallocate -l 2G /swapfile
  chmod 600 /swapfile
  mkswap /swapfile
  swapon /swapfile
  echo '/swapfile none swap sw 0 0' >> /etc/fstab
fi

mkdir -p /opt/edu-retain
if [ ! -d /opt/edu-retain/.git ]; then
  rm -rf /opt/edu-retain
  git clone https://github.com/RishiNarayan1/EduRetain.git /opt/edu-retain
fi

cd /opt/edu-retain
git fetch origin master
git reset --hard origin/master

perl -0pi -e 's/(const MONGODB_URI = .*?;\n)/$1\nif (process.env.NODE_ENV === "production") {\n    app.set("trust proxy", 1);\n}\n/s' backend/server.cjs
perl -0pi -e 's/secure: process\.env\.NODE_ENV === .production.,/secure: process.env.SESSION_COOKIE_SECURE === "true" || (\n            process.env.NODE_ENV === "production" \&\& process.env.SESSION_COOKIE_SECURE !== "false"\n        ),/' backend/server.cjs

cat > Dockerfile <<'DOCKERFILE'
FROM node:22-alpine AS build

WORKDIR /app

COPY package*.json ./
RUN npm ci

COPY . .
RUN npm run build

FROM node:22-alpine AS runtime

ENV NODE_ENV=production
ENV PORT=8080

WORKDIR /app

COPY package*.json ./
RUN npm ci --omit=dev && npm cache clean --force

COPY backend ./backend
COPY src/data ./src/data
COPY --from=build /app/dist ./dist

EXPOSE 8080

CMD ["node", "backend/server.cjs"]
DOCKERFILE

SESSION_SECRET="$(openssl rand -hex 32)"

docker network create edu-retain-net || true
docker volume create edu-retain-mongo
docker rm -f edu-retain-app edu-retain-mongo || true

docker run -d \
  --name edu-retain-mongo \
  --network edu-retain-net \
  --restart unless-stopped \
  -v edu-retain-mongo:/data/db \
  mongo:7 --wiredTigerCacheSizeGB 0.25

until docker exec edu-retain-mongo mongosh --quiet --eval "db.adminCommand('ping').ok" | grep -q 1; do
  sleep 3
done

docker build -t edu-retain-app .
docker run --rm \
  --network edu-retain-net \
  -e MONGODB_URI=mongodb://edu-retain-mongo:27017/EDU_RETAIN \
  edu-retain-app npm run migrate || true
docker run --rm \
  --network edu-retain-net \
  -e MONGODB_URI=mongodb://edu-retain-mongo:27017/EDU_RETAIN \
  edu-retain-app npm run import-csv
docker run -d \
  --name edu-retain-app \
  --network edu-retain-net \
  --restart unless-stopped \
  -p 80:8080 \
  -e NODE_ENV=production \
  -e PORT=8080 \
  -e MONGODB_URI=mongodb://edu-retain-mongo:27017/EDU_RETAIN \
  -e SESSION_SECRET="${SESSION_SECRET}" \
  -e SESSION_COOKIE_SECURE=false \
  -e CLIENT_URL=http://localhost \
  edu-retain-app
