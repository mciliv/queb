.PHONY: help setup start stop restart logs health test clean update

help: ## Show this help message
	@echo 'Usage: make [target]'
	@echo ''
	@echo 'Available targets:'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  %-15s %s\n", $$1, $$2}'

setup: ## Initial setup - copy .env.example to .env
	@if [ ! -f .env ]; then \
		cp .env.example .env; \
		echo "✅ Created .env file from .env.example"; \
		echo "⚠️  Please edit .env and add your API keys"; \
	else \
		echo "⚠️  .env already exists, skipping"; \
	fi

start: ## Start the LiteLLM proxy server
	docker compose up -d
	@echo "✅ LiteLLM proxy started at http://localhost:4000"
	@echo "📊 Grafana dashboard at http://localhost:3000"
	@echo "📈 Prometheus at http://localhost:9090"

stop: ## Stop the proxy server
	docker compose down
	@echo "✅ LiteLLM proxy stopped"

restart: ## Restart the proxy server
	docker compose restart litellm
	@echo "✅ LiteLLM proxy restarted"

logs: ## View proxy logs (Ctrl+C to exit)
	docker compose logs -f litellm

health: ## Check proxy health
	@curl -s http://localhost:4000/health | jq . || echo "❌ Proxy not responding"

models: ## List available models
	@curl -s http://localhost:4000/models \
		-H "Authorization: Bearer $$(grep LITELLM_MASTER_KEY .env | cut -d '=' -f2)" \
		| jq '.data[] | {id: .id, created: .created}' || echo "❌ Failed to fetch models"

test: ## Test the proxy with a simple request
	@echo "Testing with GPT-4o..."
	@curl -s http://localhost:4000/chat/completions \
		-H "Authorization: Bearer $$(grep LITELLM_MASTER_KEY .env | cut -d '=' -f2)" \
		-H "Content-Type: application/json" \
		-d '{"model":"gpt-4o","messages":[{"role":"user","content":"Say hello!"}]}' \
		| jq . || echo "❌ Test failed"

clean: ## Stop and remove all containers, volumes, and data
	docker compose down -v
	@echo "✅ All containers and volumes removed"

update: ## Update to latest LiteLLM image
	docker compose pull
	docker compose up -d
	@echo "✅ Updated to latest LiteLLM version"

stats: ## Show resource usage
	docker stats --no-stream litellm postgres prometheus grafana

shell: ## Open a shell in the LiteLLM container
	docker compose exec litellm /bin/sh
