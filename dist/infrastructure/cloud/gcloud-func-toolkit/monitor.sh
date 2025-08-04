#!/bin/bash

# GCloud Function Monitoring Utilities
# ===================================

# Setup function monitoring
setup_monitoring() {
    local function_name="$1"
    
    log_info "Setting up monitoring for function: $function_name"
    
    # Enable Cloud Monitoring API
    log_info "Enabling Cloud Monitoring API..."
    gcloud services enable monitoring.googleapis.com --project="$PROJECT_ID"
    
    # Create monitoring workspace if needed
    log_info "Setting up monitoring workspace..."
    gcloud monitoring workspaces create --project="$PROJECT_ID" --display-name="$PROJECT_NAME Monitoring" 2>/dev/null || true
    
    log_info "✅ Monitoring setup complete for: $function_name"
}

# Check function logs
check_logs() {
    local function_name="$1"
    local lines="${2:-50}"
    
    log_info "Checking logs for function: $function_name"
    
    gcloud functions logs read "$function_name" \
        --region="$REGION" \
        --project="$PROJECT_ID" \
        --limit="$lines" \
        --format="table(timestamp,severity,textPayload)"
}

# Setup alerting
setup_alerts() {
    local function_name="$1"
    
    log_info "Setting up alerts for function: $function_name"
    
    # Create alert policy for function errors
    local policy_name="$function_name-error-alert"
    
    log_info "Creating error alert policy..."
    gcloud alpha monitoring policies create \
        --policy-from-file=- \
        --project="$PROJECT_ID" << EOF
displayName: "$function_name Error Alert"
conditions:
- displayName: "$function_name Error Rate"
  conditionThreshold:
    filter: 'resource.type="cloud_function" AND resource.labels.function_name="$function_name" AND severity>=ERROR'
    comparison: COMPARISON_GREATER_THAN
    thresholdValue: 0
    duration: 300s
EOF
    
    log_info "✅ Alert policy created: $policy_name"
}

# Get function metrics
get_metrics() {
    local function_name="$1"
    local duration="${2:-1h}"
    
    log_info "Getting metrics for function: $function_name"
    
    # Get execution count
    log_info "Execution count (last $duration):"
    gcloud monitoring metrics list --filter="metric.type:cloudfunctions.googleapis.com/function/execution_count" --project="$PROJECT_ID"
    
    # Get execution time
    log_info "Execution time (last $duration):"
    gcloud monitoring metrics list --filter="metric.type:cloudfunctions.googleapis.com/function/execution_time" --project="$PROJECT_ID"
    
    # Get memory usage
    log_info "Memory usage (last $duration):"
    gcloud monitoring metrics list --filter="metric.type:cloudfunctions.googleapis.com/function/memory_usage" --project="$PROJECT_ID"
}

# Check function health
check_health() {
    local function_name="$1"
    
    log_info "Checking health for function: $function_name"
    
    # Check function status
    local status=$(gcloud functions describe "$function_name" --region="$REGION" --project="$PROJECT_ID" --format="value(status)" 2>/dev/null)
    
    if [ -n "$status" ]; then
        case "$status" in
            "ACTIVE")
                log_info "✅ Function is active"
                ;;
            "DEPLOYING")
                log_info "⏳ Function is deploying"
                ;;
            "FAILED")
                log_error "❌ Function deployment failed"
                return 1
                ;;
            *)
                log_warn "⚠️ Function status: $status"
                ;;
        esac
    else
        log_error "❌ Function not found"
        return 1
    fi
    
    # Check recent errors
    local error_count=$(gcloud functions logs read "$function_name" --region="$REGION" --project="$PROJECT_ID" --limit=100 --format="value(severity)" | grep -c "ERROR" || echo "0")
    
    if [ "$error_count" -gt 0 ]; then
        log_warn "⚠️ Found $error_count errors in recent logs"
    else
        log_info "✅ No recent errors found"
    fi
    
    return 0
}

# Setup logging sink
setup_logging_sink() {
    local function_name="$1"
    local sink_name="$function_name-logs"
    
    log_info "Setting up logging sink for function: $function_name"
    
    # Create BigQuery dataset for logs
    local dataset_name="${function_name}_logs"
    bq mk --dataset "$PROJECT_ID:$dataset_name" 2>/dev/null || true
    
    # Create logging sink
    log_info "Creating logging sink..."
    gcloud logging sinks create "$sink_name" \
        "bigquery.googleapis.com/projects/$PROJECT_ID/datasets/$dataset_name" \
        --log-filter="resource.type=cloud_function AND resource.labels.function_name=$function_name" \
        --project="$PROJECT_ID"
    
    log_info "✅ Logging sink created: $sink_name"
} 